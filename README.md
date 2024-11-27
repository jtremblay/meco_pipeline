## Release notes
This repository is the sucessor of the nrc_pipeline (https://bitbucket.org/jtremblay514/nrc_pipeline_public) project. It holds a suite of bioinformatics pipelines geared for microbial ecology sequence (i.e. NGS) data processing.

## Summary
The meco pipeline's (this repository) core libraries were forked from the GenPipes (https://doi.org/10.1101/459552) repository (formely the mugqic_pipeline repository) on August 16th 2014. It has since then been extensively developped and used for many projects. Notably, it holds two core pipelines:
  AmpliconTagger has been published in GigaScience - https://doi.org/10.1093/gigascience/giz146
  ShotgunMG has been published in Briefings in Bioinformatics - https://doi.org/10.1093/bib/bbac443. ShotgunMG is also available as a standalone Nextflow pipeline (https://github.com/jtremblay/ShotgunMG).

## Requirements
meco pipeline requires many third party software. More details in the user guides mentionned below.

This pipeline infrastructure relies on its companion meco_tools repository (https://github.com/jtremblay/meco_tools).

The pipeline will generate jobs that will need to have a proper installation of the following softwares using the Environment Module System (Lmod).

Software installation scripts are available in the meco_resources repository.
https://github.com/jtremblay/meco_resources

## Reference and help
We provide a Ubuntu 22.04 image that contains all that is needed to run both ShotgunMG and AmpliconTagger. It contains the third party softwares, modules, databases and a subset of sequencing libraries for testing purposes.

https://hub.docker.com/r/julio514/ubuntu

**An exhaustive user guide is available here for AmpliconTagger: [User Guide](http://jtremblay.github.io/amplicontagger_guide.html)**


**And here for ShotgunMG: [User Guide](http://jtremblay.github.io/shotgunmg_guide.html)**

## Pipeline diagrams
http://jtremblay.github.io/PipelineViewer/amplicontagger.html


http://jtremblay.github.io/PipelineViewer/shotgunmg.html

## LICENSE AND COPYRIGHT
This software is Copyright (c) 2024 Julien Tremblay and is freely available for use without any warranty under the same license as GenPipes itself. Refer to wrapped tools for their credits and license information.
This license does not grant you the right to use any trademark, service mark, tradename, or logo of the Copyright Holder.
This license includes the non-exclusive, worldwide, free-of-charge patent license to make, have made, use, offer to sell, sell, import and otherwise transfer the Package with respect to any patent claims licensable by the Copyright Holder that are necessarily infringed by the Package. If you institute patent litigation (including a cross-claim or counterclaim) against any party alleging that the Package constitutes direct or contributory patent infringement, then this Artistic License to you shall terminate on the date that such litigation is filed.
Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES. THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

