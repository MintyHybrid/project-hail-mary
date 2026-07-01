# potytools 0.1.0

* Initial release.
* GenBank file parser (`parse_genbank_file()`, `load_genbank_folder()`) with no external dependencies.
* Codon usage analysis: `calculate_rscu()`, `calculate_enc()`, `calculate_gc_content()`, `compare_codon_usage()`, `check_wobble_bias()`.
* Host detection and classification for monocot/dicot-infecting potyviruses.
* ELM short linear motif search via REST API (`elm_batch_search()`, `elm_search_api()`).
* Motif-flanking region extraction (`extract_motif_flanks()`).
* NCBI sequence fetching with custom coordinates (`fetch_custom_sequences()`).
* Concatenated fragment builder for batch motif analysis (`build_concat()`).
