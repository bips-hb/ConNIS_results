```bash
#!/bin/bash

# Datei, aus der die Zeilen extrahiert werden sollen
input_file="14028s_genebank.txt"
# Datei, in die die extrahierten Zeilen geschrieben werden sollen
output_file="genes_and_positions.txt"

# Leere die Ausgabedatei, falls sie bereits existiert
> "$output_file"

# Durchsuche die Datei nach den gewünschten Mustern
grep -E '    gene    |/gene=|/locus_tag=|/gene_synonym' "$input_file" > "$output_file"

```