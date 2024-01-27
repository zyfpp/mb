#download metadata from ncbicli
datasets summary genome taxon human

#extract data from json
jq -r '.reports[] | [.accession, (.assembly_info.biosample.attributes[] | select(.name == "collection_date").value? // "N/A")] | @csv' meta.json > year.csv #collection date
jq -r '.reports[] | [.accession, (.checkm_info.completeness // "N/A")] | @csv' meta.json > completeness_data.csv #completeness
jq -r '.reports[] | {accession: .accession, host: (.assembly_info.biosample.attributes[] | select(.name == "host").value? // "N/A")} | [.accession, .host] | @csv' meta.json > hosts.csv #host
jq -r '.reports[] | [.accession, (.assembly_info.biosample.attributes[] | select(.name == "geo_loc_name").value? // "N/A")] | @csv' your_file.json > geo.csv #geo
jq -r '.reports[] | [.accession, .average_nucleotide_identity.submitted_organism // "N/A"] | @csv' meta.json > spc.csv #organism species
jq -r '.reports[] | [.accession, .checkm_info.contamination] | @csv' meta.json > contamination.csv #contamination
