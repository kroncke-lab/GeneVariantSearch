# PubMed Variant Annotation Validator

## Overview
A Streamlit web application that validates genetic variant annotations against their source publications. Upload a CSV with annotated variants and PMIDs, and the tool will verify whether the annotations match what's actually reported in the papers.

## Features
- **CSV Upload**: Supports various CSV formats with auto-detection
- **Smart Validation**: Fetches full-text papers including supplements
- **AI-Powered Analysis**: Uses Gemini/Claude/OpenAI to compare annotations vs. reported data
- **Clear Results**: Returns CORRECT ✓ / INCORRECT ✗ / UNSURE ? with explanations
- **Export**: Download validated CSV with results

## How to Use

### 1. Prepare Your CSV
Your CSV should include:
- **PMID** column (required) - PubMed ID of the source paper
- **Variant** column(s) - HGVS notation, protein name, or other variant identifier
- **Phenotype/Diagnosis** column(s) - Clinical features, disease type, outcomes

### 2. Upload and Configure
1. Select your AI model (Gemini is free)
2. Add your API key in sidebar or Replit Secrets
3. Upload your CSV file
4. Map columns (PMID, Variant, Phenotype)

### 3. Run Validation
1. Choose how many rows to validate
2. Click "Start Validation"
3. Review results and download validated CSV

## Example CSV Formats

### Format 1: Cancer Variants
```csv
Source (URL only),Variant Name (HGVS DNA),Variant name (protein),Case #,Control,Cancer type,Odds ratio,P-value,PMID
https://...,c.123A>G,p.Tyr54Asn,45,120,Breast Cancer,2.1,0.003,19841300
```

### Format 2: Cardiac Variants
```csv
Variant,LQT,HOM,SMPT,ASD,SCD,Patient Locations,Patient ethnic background,PMID,Journal,Author,Year
p.Tyr54Asn,Yes,No,No,No,Yes,USA,Caucasian,19841300,Circulation,Smith,2009
```

## Validation Logic

For each row, the app:
1. Fetches the paper by PMID (abstract + full text + supplements)
2. Extracts what's actually reported about the variant
3. Compares annotated variant/phenotype vs. reported data
4. Returns:
   - **CORRECT**: Annotation matches paper
   - **INCORRECT**: Clear mismatch or contradiction
   - **UNSURE**: Insufficient information or ambiguous

## Requirements
- Python 3.8+
- API key for Gemini (free), Claude, or OpenAI
- Internet connection for PubMed access

## Running the App
```bash
streamlit run validator_app/app.py --server.port 5001
```

## API Keys
Get your free API key:
- **Gemini**: https://aistudio.google.com/apikey (recommended)
- **Claude**: https://console.anthropic.com
- **OpenAI**: https://platform.openai.com

Add to Replit Secrets:
- `GEMINI_API_KEY`
- `ANTHROPIC_API_KEY` 
- `OPENAI_API_KEY`

## GitHub Repository
This app can be pushed to a separate GitHub repository. To export:
1. Copy the `validator_app/` folder
2. Create a new GitHub repo
3. Push the validator_app contents

Or use Replit's GitHub integration to sync this folder.
