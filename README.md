# XOOS Secondary Analysis

## Introduction

XOOS (SB**X** **O**ptimized **O**pen **S**ource) analysis tools are a suite of bioinformatics tools designed and optimized for processing SBX data generated on the AXELIOS platform.
The tools cover the full secondary-analysis path for SBX-Duplex (SBX-D) data — from demultiplexing and adapter trimming, through read collapsing and alignment metrics, to variant calling — and can be composed into end-to-end pipelines or run as standalone tools.

The tools support several analysis applications on SBX-D WGS data:

- **Germline WGS** — small variants (SNVs and indels), copy-number events, and repeat expansions (STR).
- **Somatic Tumor/Normal (TN) WGS** — paired tumor/normal small-variant and copy-number analysis.
- **cfDNA WGS** — circulating tumor DNA reporting against a matched Tumor/Normal result.

Individual tools are at different maturity levels (some applications are in alpha or beta); see each tool's documentation in the navigation sidebar for current chemistry and application support.

As an **open source solution**, XOOS is freely available for the scientific community to adopt, customize, and enhance.
We encourage users to:

- **Adapt** the tools to your specific needs
- **Configure** parameters and workflows for your computational environment
- **Extend** functionality with additional tools and analysis modules
- **Contribute back** improvements, bug fixes, and new features to benefit the SBX community

All XOOS analysis tools source code is hosted on [GitHub](https://github.com/Roche-AXELIOS/XOOS), please see [Releases](https://github.com/Roche-AXELIOS/XOOS/releases) for latest release notes, and [Packages](https://github.com/orgs/Roche-AXELIOS/packages?repo_name=XOOS) for pre-built Docker images.

## Disclaimer

**For Research Use Only. Not for use in diagnostic procedures.**

Unless required by applicable law or agreed to in writing, XOOS analysis tools are provided on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied, including, without limitation, any warranties or conditions of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE.
You are solely responsible for determining the appropriateness of using or redistributing the Work and assume any risks associated with your use of XOOS analysis tools.
