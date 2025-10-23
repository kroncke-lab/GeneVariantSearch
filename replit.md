# PubMed Genetic Variant & Clinical Data Extractor

## Overview
A Streamlit web application that searches PubMed publications for genetic variants and extracts clinical data using AI-powered text analysis. The tool helps researchers quickly identify published genetic variants, their carrier genotypes, associated phenotypes, and clinical outcomes.

## Purpose
Extract structured genetic and clinical data from scientific literature, specifically:
- Genes and variants (e.g., KCNH2 p.Tyr54Asn)
- Carrier genotypes (heterozygous, homozygous, compound het)
- Phenotypic data (QT prolongation, arrhythmia, syncope)
- Optional attributes (age, sex, treatment, outcomes)

## Recent Changes
- **2025-10-23**: Initial implementation
  - PubMed API integration using Biopython
  - LLM-powered variant extraction (Claude 3 Haiku or GPT-3.5 Turbo)
  - Streamlit interface with search, analysis, and CSV export
  - Support for broad gene searches or specific variant queries

## User Preferences
- Prefers cost-effective/free LLM models for initial analysis
- Wants to search PubMed XML-formatted publications
- Needs structured data extraction for further research analysis

## Project Architecture

### Core Components
1. **pubmed_api.py**: PubMed/NCBI Entrez API integration
   - Searches publications by gene/variant
   - Retrieves XML-formatted article metadata
   - Parses titles, abstracts, authors, publication year

2. **llm_analyzer.py**: AI-powered data extraction
   - Supports Anthropic Claude 3 Haiku (cheapest option)
   - Supports OpenAI GPT-3.5 Turbo
   - Extracts variants, genotypes, phenotypes, clinical data
   - Returns structured JSON responses

3. **app.py**: Streamlit web interface
   - Gene/variant input form
   - Model selection (Anthropic/OpenAI)
   - Progress tracking during analysis
   - Results display with metrics
   - CSV export functionality

### Dependencies
- streamlit: Web interface
- biopython: PubMed API access
- anthropic: Claude AI integration
- openai: GPT integration
- pandas: Data structuring and export

### API Keys Required
- `ANTHROPIC_API_KEY` or `OPENAI_API_KEY` (user choice)
- PubMed requires email but no API key

### Workflow
- Command: `streamlit run app.py --server.port 5000`
- Port: 5000 (webview)

## How to Use
1. Select AI model (Anthropic or OpenAI)
2. Enter gene name (required)
3. Enter specific variant (optional - leave empty to find all variants)
4. Set maximum articles to analyze
5. Click "Search PubMed & Analyze"
6. Review extracted data in table
7. Download results as CSV

## Example Searches
- Gene: `KCNH2`, Variant: `p.Tyr54Asn` → Find specific variant
- Gene: `BRCA1`, Variant: empty → Find all BRCA1 variants
- Gene: `SCN5A`, Variant: empty → Comprehensive variant discovery
