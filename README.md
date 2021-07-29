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