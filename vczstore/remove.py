import logging
from contextlib import nullcontext

import numpy as np
from vcztools.utils import array_dims, open_zarr, search

from vczstore.utils import missing_val, variant_chunk_slices, variants_progress

logger = logging.getLogger(__name__)


def remove(vcz, sample_id, *, show_progress=False, zarr_backend_storage=None):
    """Remove a sample from vcz and overwrite with missing data"""

    if zarr_backend_storage == "icechunk":
        from vczstore.icechunk_utils import icechunk_transaction

        cm = icechunk_transaction(vcz, "main", message="remove")
    else:
        cm = nullcontext(vcz)

    with cm as vcz:
        root = open_zarr(vcz, mode="r+", zarr_backend_storage=zarr_backend_storage)
        n_variants = root["variant_contig"].shape[0]
        all_samples = root["sample_id"][:]

        # find index of sample to remove
        unknown_samples = np.setdiff1d(sample_id, all_samples)
        if len(unknown_samples) > 0:
            raise ValueError(f"unrecognised sample: {sample_id}")
        sample_selection = search(all_samples, sample_id)

        # overwrite sample data
        root["sample_id"][sample_selection] = ""
        with variants_progress(n_variants, "Remove", show_progress) as pbar:
            for v_sel in variant_chunk_slices(root):
                for var in root.keys():
                    arr = root[var]
                    if (
                        var.startswith("call_")
                        and array_dims(arr)[0] == "variants"
                        and array_dims(arr)[1] == "samples"
                    ):
                        arr[v_sel, sample_selection, ...] = missing_val(arr)
            pbar.update(v_sel.stop - v_sel.start)
