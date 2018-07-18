#!/bin/bash
${1} aln -b0 ${2} ${3} ${4} | ${1} samse ${3} - ${4}
