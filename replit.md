# PubMed Genetic Variant Tools Suite

## Overview
Two Streamlit web applications for genetic variant research:

### 1. Variant Extractor (main app)
Searches PubMed publications for genetic variants and extracts clinical data using AI-powered text analysis. Helps researchers quickly identify published genetic variants, their carrier genotypes, associated phenotypes, and clinical outcomes.

### 2. Variant Validator (validator_app)
Validates genetic variant annotations from CSV files against source PubMed publications. Compares annotated data vs. what's actually reported in papers and returns CORRECT/INCORRECT/UNSURE assessments with explanations.

## Purpose
Extract structured genetic and clinical data from scientific literature, specifically:
- Genes and variants (e.g., KCNH2 p.Tyr54Asn)
- Carrier genotypes (heterozygous, homozygous, compound het)
- Phenotypic data (QT prolongation, arrhythmia, syncope)
- Optional attributes (age, sex, treatment, outcomes)

## Apps

### Variant Extractor (`app.py`)
**Purpose:** Discover and extract variant data from PubMed literature
**Features:**
- Search by gene name and optional specific variant
- Extracts from abstracts, full-text, and supplements (PDF/Excel/Word)
- Smart filtering to find variant-relevant sections
- Supports shorthand notations (Y54N, Tyr54Asn, p.Tyr54Asn, etc.)
- CSV export of extracted data

**Run:** `streamlit run app.py --server.port 5000`

### Variant Validator (`validator_app/app.py`)
**Purpose:** Quality control - validate annotations against source papers
**Features:**
- Upload CSV with PMID, Gene, Variant, Phenotype columns
- Auto-detects columns with smart defaults
- Fetches source papers (abstracts + full text when available)
- AI compares annotated vs. reported data
- Returns CORRECT/INCORRECT/UNSURE with explanations
- CSV export with validation results

**Run:** `streamlit run validator_app/app.py --server.port 5000`

**Tip:** For best results, include a Gene column in your CSV for better full-text retrieval.

## Recent Changes
- **2025-10-27**: Added Variant Validator app
  - New app in `validator_app/` for validating CSV annotations
  - Supports optional gene column for improved full-text retrieval
  - Uses GPT-4o Mini (instead of GPT-3.5) for JSON output support
  - Auto-detects PMID/Gene/Variant/Phenotype columns
  - Returns validation status with detailed explanations
- **2025-10-24**: Added Gemini AI support
  - Added Google Gemini 2.5 Flash as AI model option (free tier available)
  - Updated UI to show Gemini as default/first choice
  - Users can now choose between Gemini, Claude, or GPT models
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
   - Supports Google Gemini 2.5 Flash (free tier available)
   - Supports Anthropic Claude 3 Haiku
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
- google-genai: Gemini AI integration
- anthropic: Claude AI integration
- openai: GPT integration
- pandas: Data structuring and export

### API Keys Required
- `GEMINI_API_KEY`, `ANTHROPIC_API_KEY`, or `OPENAI_API_KEY` (user choice)
- PubMed requires email but no API key
- Gemini offers a free tier at https://aistudio.google.com/apikey

### Workflow
- Command: `streamlit run app.py --server.port 5000`
- Port: 5000 (webview)

## How to Use
1. Select AI model (Gemini, Anthropic, or OpenAI)
2. Add your API key in Secrets (Tools → Secrets)
3. Enter gene name (required)
4. Enter specific variant (optional - leave empty to find all variants)
5. Set maximum articles to analyze
6. Click "Search PubMed & Analyze"
7. Review extracted data in table
8. Download results as CSV

## Example Searches
- Gene: `KCNH2`, Variant: `p.Tyr54Asn` → Find specific variant
- Gene: `BRCA1`, Variant: empty → Find all BRCA1 variants
- Gene: `SCN5A`, Variant: empty → Comprehensive variant discovery
