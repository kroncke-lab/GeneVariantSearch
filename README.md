# GeneVariantSearch

A PubMed literature mining tool that extracts genetic variant data and clinical information from research articles using AI analysis.

## Features

### Dual-Level Data Extraction

The tool extracts **two types** of variant data from publications:

#### 1. Individual Patient Data
- Case reports and patient narratives
- Demographics: age, sex, ethnicity/ancestry
- Genotype information (heterozygous, homozygous, compound heterozygous)
- Clinical phenotypes and diagnoses
- Treatment and outcomes

#### 2. Aggregate Study Data
- Population-level statistics from cohort studies, GWAS, and meta-analyses
- Case/control sample sizes
- Variant frequencies in cases vs controls
- Ethnicity/ancestry breakdowns (e.g., "500 EUR, 200 EAS, 100 AFR")
- Study location and institution
- Cohort names for tracking participant overlap

### Study Metadata & Deduplication

- **Study Type Detection**: Identifies case-control, cohort, GWAS, meta-analysis, family studies
- **Geographic Information**: Extracts country, city, institution, population source
- **Cohort Tracking**: Captures official cohort names (UK Biobank, CHARGE Consortium, etc.)
- **Overlap Detection**: Automatically warns when the same cohort appears in multiple studies to prevent double-counting participants

### Data Sources

- **PubMed**: Full-text search and article metadata
- **PubMed Central (PMC)**: Open-access full-text articles with inline tables
- **Supplemental Materials**: Automatically downloads and parses:
  - Excel files (.xlsx, .xls)
  - PDF documents
  - Word documents (.docx, .doc)

### AI Models Supported

- **Google Gemini 2.5 Flash** (recommended, FREE tier available)
- **Anthropic Claude 3 Haiku**
- **OpenAI GPT-3.5 Turbo**

## Installation

1. Install dependencies:
```bash
pip install -r requirements.txt
```

2. Set up API key for your chosen AI model:
```bash
# For Gemini (recommended - free tier)
export GEMINI_API_KEY="your_api_key_here"

# Or for Anthropic Claude
export ANTHROPIC_API_KEY="your_api_key_here"

# Or for OpenAI
export OPENAI_API_KEY="your_api_key_here"
```

Get a free Gemini API key at: https://aistudio.google.com/apikey

## Usage

Run the Streamlit app:
```bash
streamlit run app.py
```

### Search Options

1. **Gene Name** (required): Standard gene symbol (e.g., KCNH2, BRCA1, SCN5A)
2. **Variant** (optional): Specific variant notation (e.g., p.Tyr54Asn, c.161A>T)
3. **Target Phenotypes**: Comma-separated list of clinical features to prioritize (e.g., "QT interval, arrhythmia, sudden death")
4. **Max Articles**: Number of articles to analyze (5-200)

### Example Queries

- **Gene**: `KCNH2`, **Variant**: `p.Tyr54Asn`, **Phenotypes**: `QT interval, arrhythmia`
- **Gene**: `BRCA1`, **Phenotypes**: `breast cancer, ovarian cancer`
- **Gene**: `SCN5A` (all variants)

## Output

### Two Separate Data Tables

1. **Individual Patient Data**
   - Columns: PMID, Title, Year, Study Type, Location, Cohort, Gene, Variant, Genotype, Phenotypes, Age, Sex, Ethnicity, Treatment, Outcome, Confidence

2. **Aggregate Study Data**
   - Columns: PMID, Title, Year, Study Type, Location, Cohort, Gene, Variant, Cases w/ Variant, Total Cases, Controls w/ Variant, Total Controls, Ethnicity Breakdown, AF Cases, AF Controls, Phenotype, Notes, Confidence

### Export Options

- Download **Individual Patient Data** as CSV
- Download **Aggregate Study Data** as CSV
- Both exports timestamped for version tracking

### Deduplication Alerts

The tool automatically detects when multiple articles reference the same cohort:

```
⚠️ Potential Participant Overlap Detected
The following cohorts appear in multiple studies. Participants may be duplicated:
- UK Biobank: PMID 12345678 (2020), PMID 23456789 (2022)
- CHARGE Consortium: PMID 34567890 (2019), PMID 45678901 (2021)
```

Review flagged studies carefully to avoid counting the same individuals multiple times.

## Use Cases

### Research Applications

1. **Meta-Analysis Preparation**: Collect case/control counts from multiple studies
2. **Variant Discovery**: Find all reported carriers of a specific variant
3. **Phenotype Association**: Identify clinical presentations across populations
4. **Population Genetics**: Track variant frequencies across ethnicities
5. **Literature Review**: Systematic extraction of genetic evidence

### Avoiding Duplicate Data

When analyzing aggregate data:
- Check deduplication warnings for cohort overlap
- Note recruitment periods to identify temporal overlap
- Track geographic locations and institutions
- Review author lists for potential participant reuse in follow-up studies

## Technical Details

### Data Processing Pipeline

```
PubMed Search
    ↓
Fetch Article Metadata + Full Text (PMC)
    ↓
Download Supplemental Materials (Excel/PDF/Word)
    ↓
Intelligent Filtering (gene/variant-relevant sections)
    ↓
AI Analysis (structured JSON extraction)
    ↓
Dual Output: Individual Patients + Aggregate Statistics
    ↓
Deduplication Warnings + CSV Export
```

### Key Components

- **`app.py`**: Streamlit UI with dual-table display
- **`pubmed_api.py`**: PubMed/PMC data retrieval, supplement parsing
- **`llm_analyzer.py`**: AI-powered extraction with dual-level schema

## Limitations

- **PMC Access**: Only open-access articles have full text + supplements; others use abstracts only
- **AI Accuracy**: Results depend on LLM performance; verify critical findings manually
- **No Variant Validation**: Doesn't validate variants against reference databases (use ClinVar/gnomAD for that)
- **Deduplication**: Overlap detection based on cohort names only; manual review still required
- **Rate Limiting**: NCBI API limits apply (built-in 0.34s delays between requests)

## Future Enhancements

- Integration with ClinVar for clinical significance
- gnomAD population frequency lookup
- HGVS variant normalization
- HPO phenotype standardization
- Advanced duplicate detection using author overlap and recruitment dates

## License

See LICENSE file for details.

## Citation

If you use this tool in your research, please cite the repository and the AI models used for extraction.

## Support

For issues or feature requests, please open an issue on GitHub.
