# Plots bigwig signal given one or more gene coordinates

# Installation
Install the required packages:
- [samtools](https://www.htslib.org/)
- [pysam](https://pypi.org/project/pysam/)
- [pybigwig](https://pypi.org/project/pyBigWig/)

Then, run: ```python setup.py install```
# Usage:

```
bigwig-plotter \
--bam examples/LARP6.CTRL_IN1.umi.r1.fq.genome-mapped.sorted.bam \
--chrom_sizes examples/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes \
-r examples/GAPDH_example.bed \
--width 12 \
--height 2 \
--rpm \
--output examples/image_rpm_scaled.svg

bigwig-plotter \
--bam examples/LARP6.CTRL_IN1.umi.r1.fq.genome-mapped.sorted.bam \
--chrom_sizes examples/GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes \
-r examples/RBFOX2_example.bed \
--width 12 \
--height 1.5 \
--rpm \
--output examples/image_unscaled.svg
```