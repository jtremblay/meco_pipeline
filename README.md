## Release notes
v1.1: Updated all the Python 2.7 scripts to Python 3 (v3.9.0). Replaced legacy QC method (tagsQC.pl) from caf_tools with a combination of bbduk and pTrimmer which results in significant gains in processing speed.


## LICENSE AND COPYRIGHT
This software is Copyright (c) 2023 INRS - Centre Armand-Frappier and is freely available for use without any warranty under the same license as GenPipes itself. Refer to wrapped tools for their credits and license information.
This license does not grant you the right to use any trademark, service mark, tradename, or logo of the Copyright Holder.
This license includes the non-exclusive, worldwide, free-of-charge patent license to make, have made, use, offer to sell, sell, import and otherwise transfer the Package with respect to any patent claims licensable by the Copyright Holder that are necessarily infringed by the Package. If you institute patent litigation (including a cross-claim or counterclaim) against any party alleging that the Package constitutes direct or contributory patent infringement, then this Artistic License to you shall terminate on the date that such litigation is filed.
Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES. THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Summary
The CAF's pipeline's (this repository) core libraries were forked from the GenPipes (https://doi.org/10.1101/459552) repository (formely the mugqic_pipeline repository) on August 16th 2014.
AmpliconTagger has been published in GigaScience - https://doi.org/10.1093/gigascience/giz146
ShotgunMG will be available soon.

## Requirements
The software from this repository requires Python, Perl and R and has been tested with the following versions:
Python/3.9.0; Perl/5.26; R/3.6.0

This pipeline infrastructure relies on its companion caf_tools repository (https://bitbucket.org/jtremblay514/caf_tools/src/1.2/).

The pipeline will generate jobs that will need to have a proper installation of the following softwares using the Environment Module System (Lmod).
pynast/1.2.2; perl/5.26.0; rdp_classifier/2.5; fasttree/2.1.10; FLASH/1.2.11; qiime/1.9.1; duk/1.051; dnaclust/3; fastx/0.0.13.2; python/3.9.0;
python/3.6.5; R/3.6.0; java/jdk1.8.0_144; blast/2.6.0+; picrust/1.1.0; deblur/1.0.4; vsearch/2.7.1; bbmap/38.11

Software installation scripts are available in the caf_resources repository.
https://bitbucket.org/jtremblay514/caf_resources/src/1.1/


## Reference and help
We provide a CentOS-7 image that contains all that is needed to run AmpliconTagger. It contains the third party softwares, modules, databases and a subset of sequencing libraries for testing purposes.

https://hub.docker.com/r/julio514/centos

**An exhaustive user guide is available here: [User Guide](http://jtremblay.github.io/amplicontagger_guide.html)**


## Pipeline diagram
http://jtremblay.github.io/PipelineViewer/amplicontagger.html
