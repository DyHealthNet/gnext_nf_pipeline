"""
Create json files which can be used to render Manhattan plots.

Extracted from locuszoom-hosted: d726fb2
"""
# NOTE: `qval` means `-log10(pvalue)`

# TODO: optimize binning for fold@20 view.
#       - if we knew the max_qval before we started (eg, by running qq first), it would be very easy.
#       - at present, we set qval bin size well for the [0-40] range but not for variants above that.
# TODO: combine with QQ?

import collections
import heapq
import logging
import math

from zorp.parsers import BasicVariant

logger = logging.getLogger(__name__)


class MaxPriorityQueue:
    """
    .pop() returns the item with the largest priority.
    .popall() iteratively .pop()s until empty.
    priorities must be comparable.
    `item` can be anything.
    """
    # TODO: check if this is slower than blist-based MaxPriorityQueue, for ~500 items
    # Note: `ComparesFalse()` is used to prevent `heapq` from comparing `item`s to eachother.
    #       Even if two priorities are equal, `ComparesFalse() <= ComparesFalse()` will be `False`, so `item`s won't
    #       be compared.
    class ComparesFalse:
        __eq__ = __lt__ = __gt__ = lambda s, o: False

    def __init__(self):
        self._q = []  # a heap-property-satisfying list like [(priority, ComparesFalse(), item), ...]

    def add(self, item, priority):
        heapq.heappush(self._q, (priority, MaxPriorityQueue.ComparesFalse(), item))

    def add_and_keep_size(self, item, priority, size, popped_callback):
        """
        This is a min-heap data structure, but the "priorities" we receive are in terms of -log10p (larger is better)
        In order to satisfy the "smaller values are best" property, flip the sign on the priority as we add to heap
        """
        if len(self._q) < size:
            self.add(item, priority)
        else:
            if priority > self._q[0][0]:  # if new priority < the biggest priority in the heap, switch them
                _, _, item = heapq.heapreplace(self._q, (priority, MaxPriorityQueue.ComparesFalse(), item))
            popped_callback(item)

    def pop(self):
        _, _, item = heapq.heappop(self._q)
        return item

    def __len__(self):
        return len(self._q)

    def pop_all(self):
        while self._q:
            yield self.pop()


class Binner:
    """Manhattan plot binner class"""
    def __init__(self, *,
                 peak_neg_log_pval_threshold: float = 6.0,
                 peak_sprawl_dist: int = int(200e3),
                 peak_max_count: int = 500,
                 num_unbinned: int = 500,
                 bin_length: int = int(3e6)):
        # Instance configuration
        self._peak_neg_log_pval_threshold = peak_neg_log_pval_threshold
        self._peak_sprawl_dist = peak_sprawl_dist
        self._peak_max_count = peak_max_count
        self._num_unbinned = num_unbinned
        self._bin_length = bin_length

        # Internal storage
        self._peak_best_variant = None
        self._peak_last_chrpos = None
        self._peak_pq = MaxPriorityQueue()
        self._unbinned_variant_pq = MaxPriorityQueue()
        self._bins = collections.OrderedDict()  # like {<chrom>: {<pos // bin_length>: [{chrom, startpos, qvals}]}}
        self._qval_bin_size = 0.05  # this makes 200 bins for the minimum-allowed y-axis covering 0-10

    def process_variant(self, variant: BasicVariant):
        """
        There are 3 types of variants:
          a) If the variant starts or extends a peak and has a stronger pval than the current `peak_best_variant`:
             1) push the old `peak_best_variant` into `unbinned_variant_pq`.
             2) make the current variant the new `peak_best_variant`.
          b) If the variant ends a peak, push `peak_best_variant` into `peak_pq` and push the current variant into
            `unbinned_variant_pq`.
          c) Otherwise, just push the variant into `unbinned_variant_pq`.
        Whenever `peak_pq` exceeds the size `conf.manhattan_peak_max_count`, push its member with the weakest pval into
            `unbinned_variant_pq`.
        Whenever `unbinned_variant_pq` exceeds the size `conf.manhattan_num_unbinned`, bin its member with the
            weakest pval.
        So, at the end, we'll have `peak_pq`, `unbinned_variant_pq`, and `bins`.
        """

        # TODO: Internally, PheWeb binner relies on the data being mutable.
        #   Hence the container type defines what fields we use, but internally variants must be represented as dicts.
        #   This is low hanging fruit for optimization in the future (eg using a mutable, non-slots container
        #   in the parser by default)
        variant_dict = variant.to_dict()
        variant_dict['pvalue'] = variant.pvalue  # derived property

        if not math.isinf(variant_dict['pvalue']):
            # Determine bin size based on variants with finite values (eg not pvalue underflow)
            qval = variant_dict['neg_log_pvalue']
            if qval > 40:
                # this makes 200 bins for a y-axis extending past 40 (but folded so that the lower half is 0-20)s
                self._qval_bin_size = 0.2
            elif qval > 20:
                self._qval_bin_size = 0.1  # this makes 200-400 bins for a y-axis extending up to 20-40.

        if variant_dict['neg_log_pvalue'] > self._peak_neg_log_pval_threshold:  # part of a peak
            if self._peak_best_variant is None:  # open a new peak
                self._peak_best_variant = variant_dict
                self._peak_last_chrpos = (variant_dict['chrom'], variant_dict['pos'])
            elif self._peak_last_chrpos[0] == variant_dict['chrom'] \
                    and self._peak_last_chrpos[1] + self._peak_sprawl_dist > variant_dict['pos']:
                # If this new position is near the previous top hit, extend current peak. I *think* the spec calls for
                #   only a few top hits in a given window (not everything in a wide peak), so all the lower ones in
                #   that peak become binned
                self._peak_last_chrpos = (variant_dict['chrom'], variant_dict['pos'])
                if variant_dict['neg_log_pvalue'] <= self._peak_best_variant['neg_log_pvalue']:
                    self._maybe_bin_variant(variant_dict)
                else:
                    self._maybe_bin_variant(self._peak_best_variant)
                    self._peak_best_variant = variant_dict
            else:  # close old peak and open new peak
                self._maybe_peak_variant(self._peak_best_variant)
                self._peak_best_variant = variant_dict
                self._peak_last_chrpos = (variant_dict['chrom'], variant_dict['pos'])
        else:
            self._maybe_bin_variant(variant_dict)

    def _maybe_peak_variant(self, variant: dict):

        self._peak_pq.add_and_keep_size(variant, variant['neg_log_pvalue'],
                                        size=self._peak_max_count,
                                        popped_callback=self._maybe_bin_variant)

    def _maybe_bin_variant(self, variant):
        self._unbinned_variant_pq.add_and_keep_size(variant, variant['neg_log_pvalue'],
                                                    size=self._num_unbinned,
                                                    popped_callback=self._bin_variant)

    def _bin_variant(self, variant):
        chrom_idx = variant['chrom']  # This part differs from PheWeb
        if chrom_idx not in self._bins:
            self._bins[chrom_idx] = {}
        pos_bin_id = variant['pos'] // self._bin_length
        if pos_bin_id not in self._bins[chrom_idx]:
            self._bins[chrom_idx][pos_bin_id] = {
                'chrom': variant['chrom'],
                'startpos': pos_bin_id * self._bin_length,
                'qvals': set()
            }
        qval = self._rounded(variant['neg_log_pvalue'])
        self._bins[chrom_idx][pos_bin_id]["qvals"].add(qval)

    def get_result(self):
        self.get_result = None  # this can only be called once

        if self._peak_best_variant:
            self._maybe_peak_variant(self._peak_best_variant)

        peaks = list(self._peak_pq.pop_all())
        for peak in peaks:
            peak['peak'] = True

        unbinned_variants = list(self._unbinned_variant_pq.pop_all())
        unbinned_variants = sorted(unbinned_variants + peaks,
                                   key=(lambda variant: variant['neg_log_pvalue']),
                                   reverse=True)

        # unroll dict-of-dict-of-array `bins` into array `variant_bins`
        variant_bins = []
        for chrom_idx in sorted(self._bins.keys()):
            for pos_bin_id in sorted(self._bins[chrom_idx].keys()):
                b = self._bins[chrom_idx][pos_bin_id]
                assert len(b['qvals']) > 0
                b['qvals'], b['qval_extents'] = self._get_qvals_and_qval_extents(b['qvals'])
                b['pos'] = int(b['startpos'] + self._bin_length / 2)
                del b['startpos']
                variant_bins.append(b)

        return {
            'variant_bins': variant_bins,
            'unbinned_variants': unbinned_variants,
        }

    def _rounded(self, qval: float):
        # round down to the nearest multiple of `self._qval_bin_size`, then add 1/2 of `self._qval_bin_size` to be
        # in the middle of the bin
        x = qval // self._qval_bin_size * self._qval_bin_size + self._qval_bin_size / 2
        return round(x, 3)  # trim `0.35000000000000003` to `0.35` for convenience and network request size

    def _get_qvals_and_qval_extents(self, qvals: list):
        """The PheWeb manhattan plot UI code draws three kinds of points: significant hits (clickable),
            small (dense clusters) of hits drawn as blocky rectangles, and small (but isolated) hits drawn as single
            points.

        This code decides how to group the small hits: whether they are points or rectangles.
        """
        # Binning is something we do for UI; don't create "aggregate chunks" for infinity values (or any NaN output
        #   by rounding in a previous step)
        qvals = sorted(self._rounded(qval) for qval in qvals
                       if not math.isinf(qval) and not math.isnan(qval))
        extents = [[qvals[0], qvals[0]]]  # Initialize the extents as a single point so that we can begin grouping
        for q in qvals:
            # If this item falls near an existing bin, add to bin; else create a new bin
            if q <= (extents[-1][1] + self._qval_bin_size * 1.1):
                extents[-1][1] = q  # expand the last bin to cover this new point
            else:
                extents.append([q, q])
        rv_qvals, rv_qval_extents = [], []
        for (start, end) in extents:
            # Separate the single points from the regions
            if start == end:
                rv_qvals.append(start)
            else:
                rv_qval_extents.append([start, end])
        return (rv_qvals, rv_qval_extents)