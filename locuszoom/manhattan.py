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
from typing import Optional, Callable, Dict, List, Tuple, Any

Variant = Dict[str, Any]

BIN_LENGTH = int(3e6)

logger = logging.getLogger(__name__)

chrom_order_list = [str(i) for i in range(1, 22 + 1)] + ["X", "Y", "MT"]
chrom_order = {chrom: index for index, chrom in enumerate(chrom_order_list)}

class MaxPriorityQueue:
    # TODO: check if this is slower than blist-based MaxPriorityQueue, for ~500 items
    # Note: `ComparesFalse()` is used to prevent `heapq` from comparing `item`s to eachother.
    #       Even if two priorities are equal, `ComparesFalse() <= ComparesFalse()` will be `False`, so `item`s won't be compared.
    """
    `.pop()` returns the item with the largest priority.
    `.popall()` iteratively `.pop()`s until empty.
    priorities must be comparable.
    `item` can be anything.
    """

    class ComparesFalse:
        __eq__ = __lt__ = __gt__ = lambda s, o: False

    def __init__(self):
        self._q = []  # a heap-property-satisfying list like [(priority, ComparesFalse(), item), ...]

    def add(self, item, priority) -> None:
        heapq.heappush(self._q, (-priority, MaxPriorityQueue.ComparesFalse(), item))

    def add_and_keep_size(
        self, item, priority, size: int, popped_callback: Optional[Callable] = None
    ) -> None:
        if len(self._q) < size:
            self.add(item, priority)
        else:
            if (
                -priority > self._q[0][0]
            ):  # if the new priority isn't as big as the biggest priority in the heap, switch them
                _, _, item = heapq.heapreplace(
                    self._q, (-priority, MaxPriorityQueue.ComparesFalse(), item)
                )
            if popped_callback:
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
    def __init__(self, manhattan_num_unbinned=500, within_pheno_mask_around_peak=500_000, between_pheno_mask_around_peak=1_000_000, manhattan_peak_max_count=500, manhattan_peak_pval_threshold=1e-6, manhattan_peak_sprawl_dist=200_000, manhattan_peak_variant_counting_pval_threshold=5e-8):
        self._peak_best_variant = None
        self._peak_last_chrpos = None
        self._peak_pq = MaxPriorityQueue()
        self._unbinned_variant_pq = MaxPriorityQueue()
        self._bins = {}  # like {<chrom>: {<pos // bin_length>: [{chrom, startpos, qvals}]}}
        self._qval_bin_size = (
            0.05  # this makes 200 bins for the minimum-allowed y-axis covering 0-10
        )
        self._num_significant_in_current_peak = 0  # num variants stronger than manhattan_peak_variant_counting_pval_threshold
       
        self.manhattan_num_unbinned = manhattan_num_unbinned
        self.within_pheno_mask_around_peak = within_pheno_mask_around_peak
        self.between_pheno_mask_around_peak = between_pheno_mask_around_peak
        self.manhattan_peak_max_count = manhattan_peak_max_count
        self.manhattan_peak_pval_threshold = manhattan_peak_pval_threshold
        self.manhattan_peak_sprawl_dist = manhattan_peak_sprawl_dist
        self.manhattan_peak_variant_counting_pval_threshold = manhattan_peak_variant_counting_pval_threshold

        self.manhattan_peak_neg_log_pval_threshold = -math.log10(manhattan_peak_pval_threshold)
        self.manhattan_peak_variant_counting_neg_log_pval_threshold = -math.log10(manhattan_peak_variant_counting_pval_threshold) 

        assert (
            manhattan_peak_variant_counting_pval_threshold
            < manhattan_peak_pval_threshold
        )  # counting must be stricter than peak-extending

    def process_variant(self, variant: BasicVariant) -> None:
        """
        There are 3 types of variants:
          a) If the variant starts or extends a peak and has a stronger pval than the current `peak_best_variant`:
             1) push the old `peak_best_variant` into `unbinned_variant_pq`.
             2) make the current variant the new `peak_best_variant`.
          b) If the variant ends a peak, push `peak_best_variant` into `peak_pq` and push the current variant into `unbinned_variant_pq`.
          c) Otherwise, just push the variant into `unbinned_variant_pq`.
        Whenever `peak_pq` exceeds the size `manhattan_peak_max_count`, push its member with the weakest pval into `unbinned_variant_pq`.
        Whenever `unbinned_variant_pq` exceeds the size `manhattan_num_unbinned`, bin its member with the weakest pval.
        So, at the end, we'll have `peak_pq`, `unbinned_variant_pq`, and `bins`.
        """
        variant_dict = variant.to_dict()
        variant_dict['pval'] = variant.pvalue
        if not math.isinf(variant_dict['pval']):
            qval = variant_dict['neg_log_pvalue']
            if qval > 40:
                self._qval_bin_size = 0.2  # this makes 200 bins for a y-axis extending past 40 (but folded so that the lower half is 0-20)
            elif qval > 20:
                self._qval_bin_size = (
                    0.1  # this makes 200-400 bins for a y-axis extending up to 20-40.
                )

        if variant_dict["neg_log_pvalue"] > self.manhattan_peak_neg_log_pval_threshold:  # part of a peak
            if self._peak_best_variant is None:  # open a new peak
                self._peak_best_variant = variant_dict
                self._peak_last_chrpos = (variant_dict["chrom"], variant_dict["pos"])
                self._num_significant_in_current_peak = (
                    1
                    if variant_dict["neg_log_pvalue"]
                    < self.manhattan_peak_variant_counting_neg_log_pval_threshold
                    else 0
                )
            elif (
                self._peak_last_chrpos[0] == variant_dict["chrom"]
                and self._peak_last_chrpos[1] + self.manhattan_peak_sprawl_dist
                > variant_dict["pos"]
            ):  # extend current peak
                if (
                    variant_dict["neg_log_pvalue"]
                    < self.manhattan_peak_variant_counting_neg_log_pval_threshold
                ):
                    self._num_significant_in_current_peak += 1
                self._peak_last_chrpos = (variant_dict["chrom"], variant_dict["pos"])
                if variant_dict["neg_log_pvalue"] >= self._peak_best_variant["neg_log_pvalue"]:
                    self._maybe_bin_variant(variant_dict)
                else:
                    self._maybe_bin_variant(self._peak_best_variant)
                    self._peak_best_variant = variant_dict
            else:  # close old peak and open new peak
                self._peak_best_variant["num_significant_in_peak"] = (
                    self._num_significant_in_current_peak
                )
                self._num_significant_in_current_peak = (
                    1
                    if variant_dict["neg_log_pvalue"]
                    < self.manhattan_peak_variant_counting_neg_log_pval_threshold
                    else 0
                )
                self._maybe_peak_variant(self._peak_best_variant)
                self._peak_best_variant = variant_dict
                self._peak_last_chrpos = (variant_dict["chrom"], variant_dict["pos"])
        else:
            self._maybe_bin_variant(variant_dict)

    def _maybe_peak_variant(self, variant: Variant) -> None:
        self._peak_pq.add_and_keep_size(
            variant,
            variant["neg_log_pvalue"],
            size=self.manhattan_peak_max_count,
            popped_callback=self._maybe_bin_variant,
        )

    def _maybe_bin_variant(self, variant: Variant) -> None:
        self._unbinned_variant_pq.add_and_keep_size(
            variant,
            variant["neg_log_pvalue"],
            size=self.manhattan_num_unbinned,
            popped_callback=self._bin_variant,
        )

    def _bin_variant(self, variant: Variant) -> None:
        chrom_idx = chrom_order[variant["chrom"]]
        if chrom_idx not in self._bins:
            self._bins[chrom_idx] = {}
        pos_bin_id = variant["pos"] // BIN_LENGTH
        if pos_bin_id not in self._bins[chrom_idx]:
            self._bins[chrom_idx][pos_bin_id] = {
                "chrom": variant["chrom"],
                "startpos": pos_bin_id * BIN_LENGTH,
                "qvals": set(),
            }
        qval = self._rounded(variant['neg_log_pvalue'])
        self._bins[chrom_idx][pos_bin_id]["qvals"].add(qval)

    def get_result(self) -> Dict[str, List[Variant]]:
        # this can only be called once
        if getattr(self, "already_got_result", None):
            raise Exception()
        self.already_got_result = True

        if self._peak_best_variant:
            self._peak_best_variant["num_significant_in_peak"] = (
                self._num_significant_in_current_peak
            )
            self._maybe_peak_variant(self._peak_best_variant)

        peaks = list(self._peak_pq.pop_all())
        for peak in peaks:
            peak["peak"] = True
        unbinned_variants = list(self._unbinned_variant_pq.pop_all())
        unbinned_variants = sorted(
            unbinned_variants + peaks, key=(lambda variant: variant["neg_log_pvalue"])
        )

        # unroll dict-of-dict-of-array `bins` into array `variant_bins`
        variant_bins = []
        for chrom_idx in sorted(self._bins.keys()):
            for pos_bin_id in sorted(self._bins[chrom_idx].keys()):
                b = self._bins[chrom_idx][pos_bin_id]
                assert len(b["qvals"]) > 0
                b["qvals"], b["qval_extents"] = self._get_qvals_and_qval_extents(
                    b["qvals"]
                )
                b["pos"] = int(b["startpos"] + BIN_LENGTH / 2)
                del b["startpos"]
                variant_bins.append(b)

        return {
            "variant_bins": variant_bins,
            "unbinned_variants": unbinned_variants,
        }

    def _rounded(self, qval: float) -> float:
        # round down to the nearest multiple of `self._qval_bin_size`, then add 1/2 of `self._qval_bin_size` to be in the middle of the bin
        x = qval // self._qval_bin_size * self._qval_bin_size + self._qval_bin_size / 2
        return round(
            x, 3
        )  # trim `0.35000000000000003` to `0.35` for convenience and network request size

    def _get_qvals_and_qval_extents(
        self, qvals: List[float]
    ) -> Tuple[List[float], List[Tuple[float, float]]]:
        qvals = sorted(self._rounded(qval) for qval in qvals)
        extents = [(qvals[0], qvals[0])]
        for q in qvals:
            if q <= extents[-1][1] + self._qval_bin_size * 1.1:
                extents[-1] = (extents[-1][0], q)
            else:
                extents.append((q, q))
        rv_qvals, rv_qval_extents = [], []
        for start, end in extents:
            if start == end:
                rv_qvals.append(start)
            else:
                rv_qval_extents.append((start, end))
        return (rv_qvals, rv_qval_extents)