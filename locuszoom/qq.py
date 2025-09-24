"""
Create json files which can be used to render QQ plots.

Extracted from locuszoom-hosted: d726fb2
"""

# TODO: make gc_lambda for maf strata, and show them if they're >1.1?
# TODO: copy some changes from <https://github.com/statgen/encore/blob/master/plot-epacts-output/make_qq_json.py>

# Peter has included some original notes on the processing requirements, as follows::
# TODO: reduce QQ memory using Counter(v.qval for v in variants).
#      - but we still need to split into 4 strata using MAF. Can that be done efficiently?
#          a) we could keep balanced lists for the 4 strata, but we can only be confidently start processing variants
#          once we've read 3/4 of all variants
#          b) we could assume that, since we're sorted by chr-pos-ref-alt, MAF should be pretty randomly ordered.
#               - then we could start processing variants after reading only 10% of all variants
#               - if we're wrong, `raise StrataGuessingFailed()` and try again with sorting.
#          c) we could run manhattan before this, and make it track Counter(rounded(v.maf,2) for v in variants).

# NOTE: `qval` means `-log10(pvalue)`

import collections
import logging
import math
import typing as ty

import boltons.iterutils
import boltons.mathutils
import scipy.stats
from zorp.parsers import BasicVariant



NUM_BINS = 400
NUM_MAF_RANGES = 4
MAF_SIGFIGS = 2


logger = logging.getLogger(__name__)


Variant = collections.namedtuple('Variant', ['qval', 'maf'])


def augment_variants(variants: ty.Iterator[BasicVariant], num_samples=None):
    for var in variants:
        if var.pvalue == 0:
            # FIXME: Why does QQ plot require this stub value?
            qval = 1000  # TODO(pjvh): make an option "convert_pval0_to = [num|None]"
        else:
            qval = var.neg_log_pvalue

        af = var.maf
        if af is not None:
            af = round(af, MAF_SIGFIGS)
        yield Variant(qval=qval, maf=af)


def round_sig(x, digits):
    if x == 0:
        return 0
    elif abs(x) == math.inf or math.isnan(x):
        raise ValueError("Cannot round infinity or NaN")
    else:
        log = math.log10(abs(x))
        digits_above_zero = int(math.floor(log))
        return round(x, digits - 1 - digits_above_zero)


# TODO: Move these to unit tests
assert round_sig(0.00123, 2) == 0.0012
assert round_sig(1.59e-10, 2) == 1.6e-10


def approx_equal(a, b, tolerance=1e-4):
    return abs(a - b) <= max(abs(a), abs(b)) * tolerance


# TODO: Move these to unit tests
assert approx_equal(42, 42.0000001)
assert not approx_equal(42, 42.01)


def make_qq_stratified(variants):
    # Some variants may be missing MAF. Sort those at the end of the list (eg, lump with the common variants)
    variants = sorted(variants, key=lambda v: (v.maf is None, v.maf))

    def make_strata(idx):
        # Note: slice_indices[1] is the same as slice_indices[0] of the next slice.
        # But that's not a problem, because range() ignores the last index.
        slice_indices = (len(variants) * idx // NUM_MAF_RANGES,
                         len(variants) * (idx + 1) // NUM_MAF_RANGES)
        qvals = sorted((variants[i].qval for i in range(*slice_indices)), reverse=True)
        return {
            'maf_range': (variants[slice_indices[0]].maf,
                          variants[slice_indices[1] - 1].maf),
            'count': len(qvals),
            'qq': compute_qq(qvals),
        }

    return [make_strata(i) for i in range(NUM_MAF_RANGES)]


def make_qq_unstratified(variants, include_qq):
    qvals = sorted((v.qval for v in variants), reverse=True)
    rv = {}
    if include_qq:
        rv['qq'] = compute_qq(qvals)
    rv['count'] = len(qvals)
    rv['gc_lambda'] = {}
    for perc in ['0.5', '0.1', '0.01', '0.001']:
        gc = gc_value_from_list(qvals, float(perc))
        if math.isnan(gc) or abs(gc) == math.inf:
            logger.warning('WARNING: got gc_value {!r}'.format(gc))
        else:
            rv['gc_lambda'][perc] = round_sig(gc, 5)
    return rv


def compute_qq(qvals):
    # qvals must be in decreasing order.
    assert all(a >= b for a, b in boltons.iterutils.pairwise(qvals))

    if len(qvals) == 0:
        return []

    if qvals[0] == 0:
        logger.warning('WARNING: All pvalues are 1! How is that supposed to make a QQ plot?')
        return []

    max_exp_qval = -math.log10(0.5 / len(qvals))
    # Our QQ plot will only show `obs_qval` up to `ceil(2*max_exp_pval)`.
    # So we can drop any obs_qval above that, to save space and make sure the visible range gets all the NUM_BINS.

    # this calculation must avoid dropping points that would be shown by the calculation done in javascript.
    # `max_obs_qval` means the largest observed -log10(pvalue) that will be shown in the plot. It's usually NOT the
    # largest in the data.
    max_obs_qval = boltons.mathutils.clamp(qvals[0],
                                           lower=max_exp_qval,
                                           upper=math.ceil(2 * max_exp_qval))
    if qvals[0] > max_obs_qval:
        for qval in qvals:
            if qval <= max_obs_qval:
                max_obs_qval = qval
                break

    occupied_bins = set()
    for i, obs_qval in enumerate(qvals):
        if obs_qval > max_obs_qval:
            continue
        exp_qval = -math.log10((i + 0.5) / len(qvals))
        exp_bin = int(exp_qval / max_exp_qval * NUM_BINS)
        # TODO(pjvh): it'd be great if the `obs_bin`s started right at the lowest qval in that `exp_bin`.
        #       that way we could have fewer bins but still get a nice straight diagonal line without that
        #       stair-stepping appearance.
        obs_bin = int(obs_qval / max_obs_qval * NUM_BINS)
        occupied_bins.add((exp_bin, obs_bin))

    bins = []
    for exp_bin, obs_bin in occupied_bins:
        assert 0 <= exp_bin <= NUM_BINS, exp_bin
        assert 0 <= obs_bin <= NUM_BINS, obs_bin
        bins.append((
            exp_bin / NUM_BINS * max_exp_qval,
            obs_bin / NUM_BINS * max_obs_qval
        ))
    return {
        'bins': sorted(bins),
        'max_exp_qval': max_exp_qval,
    }


def gc_value_from_list(qvals, quantile=0.5):
    # qvals must be in decreasing order.
    assert all(a >= b for a, b in boltons.iterutils.pairwise(qvals))
    qval = qvals[int(len(qvals) * quantile)]
    pval = 10 ** -qval
    return gc_value(pval, quantile)


def gc_value(pval, quantile=0.5):
    # This should be equivalent to R: `qchisq(median_pval, df=1, lower.tail=F) / qchisq(quantile, df=1, lower.tail=F)`
    return scipy.stats.chi2.ppf(1 - pval, 1) / scipy.stats.chi2.ppf(1 - quantile, 1)


# TODO: These should be moved to unit tests
assert approx_equal(gc_value(0.49), 1.047457)  # computed using R code.
assert approx_equal(gc_value(0.5), 1)
assert approx_equal(gc_value(0.50001), 0.9999533)
assert approx_equal(gc_value(0.6123), 0.5645607)


def get_confidence_intervals(num_variants, confidence=0.95):
    one_sided_doubt = (1 - confidence) / 2

    # `variant_counts` are the numbers of variants at which we'll calculate the confidence intervals
    # any `1 <= variant_count <= num_variants-1` could be used, but scale in powers of 2 to make the CI visually smooth
    variant_counts = []
    for x in range(0, int(math.ceil(math.log2(num_variants)))):
        variant_counts.append(2 ** x)
    variant_counts.append(num_variants - 1)
    variant_counts.reverse()

    for variant_count in variant_counts:
        rv = scipy.stats.beta(variant_count, num_variants - variant_count)
        yield {
            'x': round(-math.log10((variant_count - 0.5) / num_variants), 2),
            'y_min': round(-math.log10(rv.ppf(1 - one_sided_doubt)), 2),
            'y_max': round(-math.log10(rv.ppf(one_sided_doubt)), 2),
        }