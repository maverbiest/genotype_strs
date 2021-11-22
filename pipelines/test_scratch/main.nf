#!/usr/bin/env nextflow

process testScratch {

    scratch true

    output:
    stdout result

    script:
    """
    echo "Hello World! From \$(pwd)"
    echo "\$PWD"
    """
}

result.view{ println it.trim() }