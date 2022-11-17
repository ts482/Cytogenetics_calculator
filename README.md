# Cytogenetics_calculator

# API
Ping
```
curl -i -H "Content-Type: application/json" -X POST -d '{"id":"1"}' http://localhost:5000/karyotype/ping
```

Extract
```
curl -i -H "Content-Type: application/json" -X POST -d '{"karyotype_string":"  45,X,-Y[17]/46,XY[3]   "}' http://localhost:5000/karyotype/extract
```

Send a request with FISH results
```
curl -i -H "Content-Type: application/json" -X POST -d '{"karyotype_string":"  45,X,-Y[17]/46,XY[3]   ", "fish":{"FISH_RUNX1-RUNX1T1":false}}' http://localhost:5000/karyotype/extract
```

# CLI
The CLI can be used for testing or processing full files without starting the streamlit app. By default the output file name is 
saved as [input_filename]_extracted.csv.
```
python cli.py --file=ELN2022_examples_v2.csv --version=ELN2022
```

For full list of options use
```
python cli.py --help
```


## configs
To modify the list of variants that is specifically flagged in the output, run the list of abnormalities through "setup" then add the results to the "configs" dict. E.g.
```
configs["BJH2021"] = cyto.setup(cyto.base_extraction())
```
Where cyto.base_extraction returns a list of strings. 

### Use a different config from API
Set the "version" in the API request. Compare results for different configs:

```
curl -i -H "Content-Type: application/json" -X POST -d '{"karyotype_string":"  46,xx,t(8;16)(p11;p13)[20]   ", "version": "ELN2022"}' http://localhost:5000/karyotype/extract

curl -i -H "Content-Type: application/json" -X POST -d '{"karyotype_string":"  46,xx,t(8;16)(p11;p13)[20]   ", "version": "BJH2021"}' http://localhost:5000/karyotype/extract
```

Test whether t(8;16)(p11.2;p13.3) is counted as t(8;16)(p11;p13)
```
curl -i -H "Content-Type: application/json" -X POST -d '{"karyotype_string":"  46,xx,t(8;16)(p11.2;p13.3)[20]   ", "version": "ELN2022"}' http://localhost:5000/karyotype/extract
```