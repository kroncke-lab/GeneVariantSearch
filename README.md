# PubMed Genetic Variant Tools Suite

A comprehensive toolkit for mining and validating genetic variant data from PubMed literature using AI analysis. Includes two powerful applications for variant research.

---

## 🧬 Applications

### 1. Variant Extractor (Main App)
**Purpose:** Discover and extract genetic variant data from PubMed publications

**Key Features:**
- Search PubMed by gene name and optional specific variant
- Extracts from abstracts, full-text (PMC), and supplemental materials (PDF/Excel/Word)
- Dual-level data extraction: individual patient data + aggregate study statistics
- Smart content filtering for variant-relevant sections
- Supports multiple variant notations (p.Tyr54Asn, Y54N, c.123A>G, rs12345)
- Cohort overlap detection to prevent duplicate participant counting
- CSV export for both individual and aggregate data

**Run:** 
```bash
streamlit run app.py --server.port 5000
```

### 2. Variant Validator (Quality Control App)
**Purpose:** Validate genetic variant annotations against source publications

**Key Features:**
- Upload CSV files with variant annotations
- Auto-detects PMID, Gene, Variant, and Phenotype columns
- Fetches source papers from PubMed (including full text when available)
- AI-powered comparison: annotated data vs. what's actually reported
- Returns validation status: **CORRECT** ✓ / **INCORRECT** ✗ / **UNSURE** ?
- Detailed explanations for each validation decision
- CSV export with validation results

**Run:**
```bash
streamlit run validator_app/app.py --server.port 5000
```

---

## 🚀 Quick Start

### Installation

1. **Clone this repository:**
```bash
git clone <your-repo-url>
cd <repo-name>
```

2. **Install dependencies:**
```bash
pip install streamlit biopython google-genai anthropic openai pandas pdfplumber openpyxl python-docx
```

Or use the requirements file if available:
```bash
pip install -r requirements.txt
```

3. **Set up API key** (choose one):
```bash
# Google Gemini (recommended - FREE tier available)
export GEMINI_API_KEY="your_api_key_here"

# Or Anthropic Claude
export ANTHROPIC_API_KEY="your_api_key_here"

# Or OpenAI
export OPENAI_API_KEY="your_api_key_here"
```

Get a free Gemini API key at: https://aistudio.google.com/apikey

---

## 📖 Detailed Usage

### Variant Extractor

#### Search Interface
1. **Gene Name** (required): Standard gene symbol (e.g., KCNH2, BRCA1, SCN5A)
2. **Variant** (optional): Specific variant (e.g., p.Tyr54Asn, c.161A>T)
3. **Target Phenotypes** (optional): Prioritize specific clinical features
4. **Max Articles**: Number of articles to analyze (5-200)

#### Example Queries
- Find specific variant: Gene=`KCNH2`, Variant=`p.Tyr54Asn`, Phenotypes=`QT interval, arrhythmia`
- Find all variants in gene: Gene=`BRCA1`, Phenotypes=`breast cancer`
- Broad search: Gene=`SCN5A` (leave variant blank)

#### Output Data

**Two Separate CSV Files:**

1. **Individual Patient Data**
   - Case reports and patient narratives
   - Columns: PMID, Title, Year, Study Type, Location, Gene, Variant, Genotype, Phenotypes, Age, Sex, Ethnicity, Treatment, Outcome, Confidence

2. **Aggregate Study Data**
   - Population-level statistics from cohorts and GWAS
   - Columns: PMID, Title, Year, Study Type, Location, Cohort, Gene, Variant, Cases w/ Variant, Total Cases, Controls w/ Variant, Total Controls, Ethnicity Breakdown, AF Cases, AF Controls, Confidence

**Deduplication Warnings:**
The tool automatically flags when the same cohort appears in multiple studies to prevent double-counting participants.

### Variant Validator

#### CSV Format Requirements

Your CSV should contain these columns (names will be auto-detected):

| Column Type | Examples | Required? |
|------------|----------|-----------|
| PMID | `PMID`, `PubMed ID`, `pmid` | ✅ Yes |
| Gene | `Gene`, `gene`, `Gene Symbol` | ⚠️ Optional (but recommended) |
| Variant | `Variant`, `Mutation`, `HGVS`, `Variant name (protein)` | ✅ Yes |
| Phenotype | `Phenotype`, `Diagnosis`, `Cancer type`, `Disease` | ✅ Yes |

**Example CSV Format 1 (Cancer Study):**
```csv
PMID,Gene,Variant Name (HGVS DNA),Variant name (protein),Cancer type
19841300,KCNH2,c.162A>G,p.Tyr54Asn,Breast Cancer
23456789,BRCA1,c.68_69del,p.Glu23fs,Ovarian Cancer
```

**Example CSV Format 2 (Cardiac Study):**
```csv
PMID,Variant,LQT,SCD,Patient ethnic background
19841300,KCNH2 p.Tyr54Asn,Yes,Yes,Caucasian
23456789,SCN5A p.Arg123Gln,Yes,No,Asian
```

**Tip:** Including a Gene column improves full-text retrieval and validation accuracy.

#### Validation Process

1. Upload your CSV file
2. Map columns (auto-detected with smart defaults)
3. Select number of rows to validate
4. Click "Start Validation"
5. Review results showing CORRECT/INCORRECT/UNSURE status
6. Download validated CSV with explanations

#### Validation Statuses

- **CORRECT** ✓: Both variant and phenotype match the paper (exact or semantically equivalent)
- **INCORRECT** ✗: Clear mismatch between annotation and paper
- **UNSURE** ?: Insufficient information, ambiguous reporting, or unclear match

---

## 🛠️ Technical Architecture

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
Output: CSV files with extracted/validated data
```

### Key Components

**Variant Extractor:**
- `app.py`: Streamlit UI with dual-table display
- `pubmed_api.py`: PubMed/PMC data retrieval, supplement parsing
- `llm_analyzer.py`: AI-powered extraction with dual-level schema

**Variant Validator:**
- `validator_app/app.py`: Streamlit UI for CSV upload and validation
- `validator_app/validator.py`: Validation logic and LLM integration
- Reuses `pubmed_api.py` for paper fetching

### AI Models Supported

- **Google Gemini 2.5 Flash** (recommended, FREE tier available)
- **Anthropic Claude 3 Haiku**
- **OpenAI GPT-4o Mini**

### Smart Filtering

The tool pre-filters content before AI analysis to reduce costs and improve accuracy:

**Variant Pattern Detection:**
- Standard HGVS: `p.Tyr54Asn`, `c.123A>G`, `g.456C>T`
- Three-letter shorthand: `Tyr54Asn`, `Arg123Gln`
- Single-letter shorthand: `Y54N`, `R123Q`
- rsIDs: `rs12345678`

**Table Processing:**
- Line-by-line filtering for Excel/CSV supplements
- Preserves rows containing gene names, variant patterns, or clinical keywords
- Fallback: keeps first 200 rows if no matches found

---

## 📊 Use Cases

### Research Applications

1. **Meta-Analysis Preparation**: Collect case/control counts from multiple studies
2. **Variant Discovery**: Find all reported carriers of a specific variant
3. **Phenotype Association**: Identify clinical presentations across populations
4. **Quality Control**: Validate existing variant databases against source papers
5. **Literature Review**: Systematic extraction of genetic evidence
6. **Annotation Curation**: Verify variant-phenotype associations

### Avoiding Duplicate Data

When analyzing aggregate data:
- Check deduplication warnings for cohort overlap
- Note recruitment periods to identify temporal overlap
- Track geographic locations and institutions
- Review author lists for potential participant reuse

---

## ⚠️ Limitations

### General Limitations
- **PMC Access**: Only open-access articles have full text + supplements; others use abstracts only
- **AI Accuracy**: Results depend on LLM performance; verify critical findings manually
- **Rate Limiting**: NCBI API limits apply (built-in 0.34s delays between requests)
- **No Variant Validation**: Doesn't validate variants against reference databases (use ClinVar/gnomAD for that)

### Validator-Specific Limitations
- **Gene Extraction**: Works best when CSV includes a Gene column; otherwise relies on parsing from variant strings
- **Abstract-Only Papers**: Papers without PMC full-text may have less context for validation
- **Ambiguous Cases**: Complex or contradictory papers may result in UNSURE status
- **Notation Variations**: Minor notation differences (p. vs protein name) are tolerated, but extreme variations may cause issues

---

## 🔬 Example Workflows

### Workflow 1: Discover Variants in a Gene
```bash
# Run the extractor
streamlit run app.py --server.port 5000

# Search: Gene=BRCA1, leave variant blank
# Export both individual and aggregate CSVs
# Review deduplication warnings
```

### Workflow 2: Validate Existing Annotations
```bash
# Run the validator
streamlit run validator_app/app.py --server.port 5000

# Upload your CSV with annotations
# Map columns (PMID, Gene, Variant, Phenotype)
# Validate up to 100 rows
# Download results with validation status
```

### Workflow 3: Combined Discovery + Validation
```bash
# 1. Use extractor to find variants
# 2. Export CSV results
# 3. Switch to validator
# 4. Upload the extracted CSV to verify against original papers
```

---

## 🚀 Deployment

### Running Locally
```bash
# Extractor
streamlit run app.py --server.port 5000

# Validator
streamlit run validator_app/app.py --server.port 5000
```

### Running on Replit
1. Fork this project on Replit
2. Set API key in Secrets (Tools → Secrets)
3. Run workflow: Click the Run button
4. Switch apps by changing the workflow command

### Publishing to Web
Use Replit's built-in deployment or deploy via:
- Streamlit Cloud
- Heroku
- Google Cloud Platform
- AWS

---

## 📝 File Structure

```
.
├── app.py                      # Main variant extractor app
├── pubmed_api.py              # PubMed/PMC API integration
├── llm_analyzer.py            # AI-powered data extraction
├── validator_app/             # Variant validator app
│   ├── app.py                 # Validator Streamlit UI
│   ├── validator.py           # Validation logic
│   ├── README.md              # Validator-specific docs
│   └── .streamlit/
│       └── config.toml        # Streamlit config
├── .streamlit/
│   └── config.toml            # Streamlit config
├── replit.md                  # Replit project documentation
└── README.md                  # This file
```

---

## 🆘 Troubleshooting

### Common Issues

**Problem:** "No full text available"
- **Solution:** Paper is not in PMC open access; only abstract is available. Try searching for more recent papers or check PMC directly.

**Problem:** Validator returns many UNSURE results
- **Solution:** Add a Gene column to your CSV for better full-text retrieval. Ensure PMIDs are correct and papers are accessible.

**Problem:** Excel supplements not parsing
- **Solution:** Check that openpyxl is installed: `pip install openpyxl`

**Problem:** OpenAI JSON format error
- **Solution:** Updated code uses GPT-4o Mini which supports JSON. Ensure you're using the latest version.

**Problem:** Rate limiting errors
- **Solution:** The code includes 0.34s delays. If still seeing errors, reduce max_articles or wait before retrying.

---

## 🎯 Future Enhancements

- [ ] Integration with ClinVar for clinical significance
- [ ] gnomAD population frequency lookup
- [ ] HGVS variant normalization
- [ ] HPO phenotype standardization
- [ ] Advanced duplicate detection using author overlap
- [ ] Batch processing for large validation datasets
- [ ] API endpoint for programmatic access
- [ ] Support for additional variant databases (COSMIC, LOVD)

---

## 📜 License

See LICENSE file for details.

---

## 🙏 Citation

If you use this tool in your research, please cite:
- This repository
- The AI model used for extraction (Gemini, Claude, or GPT)
- PubMed/PMC as the data source

---

## 💬 Support

- **Issues**: Open an issue on GitHub
- **Questions**: Check existing issues or start a discussion
- **Contributing**: Pull requests welcome!

---

## 🌟 Acknowledgments

Built with:
- [Streamlit](https://streamlit.io/) - Web interface
- [Biopython](https://biopython.org/) - PubMed API access
- [Google Gemini](https://ai.google.dev/) - AI extraction
- [Anthropic Claude](https://www.anthropic.com/) - AI extraction
- [OpenAI GPT](https://openai.com/) - AI extraction
- [pdfplumber](https://github.com/jsvine/pdfplumber) - PDF parsing
- [openpyxl](https://openpyxl.readthedocs.io/) - Excel parsing
- [python-docx](https://python-docx.readthedocs.io/) - Word parsing

---

**Happy variant mining! 🧬🔬**
