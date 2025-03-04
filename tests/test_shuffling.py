from bioutils.normalize import normalize


def test_shuffling():
    ref = "TAAAAAAAT"
    alt = "TATAAAAAAT"

    chrom_seq = ref
    start = 0
    end = len(ref)
    shuffle_direction = "EXPAND"

    shuffled_interval, shuffled_alleles = normalize(
        chrom_seq, interval=(start, end), alleles=(None, alt), mode=shuffle_direction
    )

    assert shuffled_interval == (2, 2)
    assert shuffled_alleles[0] == ""
    assert shuffled_alleles[1] == "T"
