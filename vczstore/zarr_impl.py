import logging

import numpy as np
import zarr
from vcztools.retrieval import variant_iter
from vcztools.utils import search
from vcztools.vcf_writer import dims

from vczstore.utils import missing_val

logger = logging.getLogger(__name__)


def append(vcz1, vcz2, *, allow_variant_subset=False, variants_mask=None):
    """Append vcz2 to vcz1 in place"""
    root1 = zarr.open(vcz1, mode="r+")
    root2 = zarr.open(vcz2, mode="r")

    # check preconditions
    if not allow_variant_subset:
        n_variants1 = root1["variant_contig"].shape[0]
        n_variants2 = root2["variant_contig"].shape[0]
        if n_variants1 != n_variants2:
            raise ValueError(
                "Stores being appended must have same number of variants. "
                f"First has {n_variants1}, second has {n_variants2}"
            )
        for field in (
            "contig_id",
            "variant_contig",
            "variant_position",
            "variant_allele",
        ):
            values1 = root1[field][:]
            values2 = root2[field][:]
            if np.any(values1 != values2):
                raise ValueError(
                    f"Stores being appended must have same values for field '{field}'"
                )
    # TODO: different check if allow_variant_subset

    # append samples
    sample_id1 = root1["sample_id"]
    sample_id2 = root2["sample_id"]

    old_num_samples = sample_id1.shape[0]
    new_num_samples = old_num_samples + sample_id2.shape[0]
    new_shape = (new_num_samples,)
    sample_id1.resize(new_shape)
    sample_id1[old_num_samples:new_num_samples] = sample_id2[:]

    if allow_variant_subset and variants_mask is None:
        variants_mask = index_variants(vcz1, vcz2)
        # TODO: store mask

    if variants_mask is None:
        variants_sel = slice(None)
    else:
        variants_sel = ~variants_mask

    # append genotype fields
    for var in root1.keys():
        if var.startswith("call_"):
            arr = root1[var]
            if arr.ndim == 2:
                new_shape = (arr.shape[0], new_num_samples)
                arr.resize(new_shape)
                arr[variants_sel, old_num_samples:new_num_samples] = root2[var][:]
            elif arr.ndim == 3:
                new_shape = (arr.shape[0], new_num_samples, arr.shape[2])
                arr.resize(new_shape)
                arr[variants_sel, old_num_samples:new_num_samples, :] = root2[var][:]
            else:
                raise ValueError("unsupported number of dims")


def remove(vcz, sample_id):
    """Remove a sample from vcz and overwrite with missing data"""
    root = zarr.open(vcz, mode="r+")
    all_samples = root["sample_id"][:]

    # find index of sample to remove
    unknown_samples = np.setdiff1d(sample_id, all_samples)
    if len(unknown_samples) > 0:
        raise ValueError(f"unrecognised sample: {sample_id}")
    selection = search(all_samples, sample_id)

    # overwrite sample data
    root["sample_id"][selection] = ""
    for var in root.keys():
        arr = root[var]
        if (
            var.startswith("call_")
            and dims(arr)[0] == "variants"
            and dims(arr)[1] == "samples"
        ):
            root[var][:, selection, ...] = missing_val(arr)


def index_variants(vcz1, vcz2):
    """Return a mask for variants of vcz2 not in vcz1"""
    root1 = zarr.open(vcz1, mode="r")
    n_variants = root1["variant_contig"].shape[0]

    # TODO: check contig IDs are identical too

    fields = ["variant_contig", "variant_position", "variant_allele"]
    it1 = variant_iter(vcz1, fields=fields)
    it2 = variant_iter(vcz2, fields=fields)
    mask = np.zeros(n_variants, dtype=bool)
    prev = None
    for i, variant in enumerate(it1):
        if prev is not None:
            if variant_is_not_after(variant, prev):
                v = prev
            else:
                # TODO: write test for this case (and change message)
                raise ValueError(f"Variant not found in VARIANTS_VCF_FILE: {prev}")
        else:
            try:
                v = next(it2)
            except StopIteration:
                v = None
        if v is not None and variant_alleles_are_equivalent(variant, v):
            # mask is False
            prev = None
        else:
            mask[i] = True
            prev = v
    return mask


def variant_is_not_after(variant1, variant2):
    """Test if variant 1 is not after variant 2 along the genome"""
    return (
        variant1["variant_contig"] == variant2["variant_contig"]
        and variant1["variant_position"] <= variant2["variant_position"]
    )


def variant_alleles_are_equivalent(variant1, variant2):
    """Test if two variants represent equivalent alleles"""
    return (
        variant1["variant_contig"] == variant2["variant_contig"]
        and variant1["variant_position"] == variant2["variant_position"]
        and np.all(variant1["variant_allele"] == variant2["variant_allele"])
    )
