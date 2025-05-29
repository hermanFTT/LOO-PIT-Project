rule compile:
    params:
         nsims=200,
         nbins=15
    output:
        output_file="results/results-{style}-{dgp}-df_{df}-beta_{beta}-n_{n}.rds",
        rates_file="results/rates-{style}-{dgp}-df_{df}-beta_{beta}-n_{n}.rds"
    script:
        "main.R"
