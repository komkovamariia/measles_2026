# Antigen reference

The antigen annotation step downloads the official [IEDB TCRMatch table](https://github.com/IEDB/TCRMatch/blob/master/data/IEDB_data.tsv) on the first run and stores it at:

```text
resources/IEDB_data.tsv
```

The file is ignored by Git because it is an external database resource. Set `MEASLES_IEDB_FILE` to use a fixed local copy for a frozen analysis release.

The primary project analysis uses exact CDR3-beta amino-acid matches only.
