# Local Installation Guide

This guide will help you set up and run the PubMed Genetic Variant Tools on your local machine.

## Quick Installation

Run this single command to install all dependencies:

```bash
pip install streamlit biopython google-genai anthropic openai pandas pdfplumber openpyxl python-docx
```

## Individual Package Installation

If you prefer to install packages one at a time:

```bash
# Core framework
pip install streamlit

# PubMed API access
pip install biopython

# AI models (install at least one)
pip install google-genai      # For Google Gemini (FREE tier available)
pip install anthropic          # For Anthropic Claude
pip install openai             # For OpenAI GPT

# Data processing
pip install pandas

# Document parsing
pip install pdfplumber         # PDF support
pip install openpyxl           # Excel support
pip install python-docx        # Word document support
```

## Setting Up API Keys

### Option 1: Environment Variables (Recommended)

**On macOS/Linux:**
```bash
# Add to ~/.bashrc or ~/.zshrc
export GEMINI_API_KEY="your_api_key_here"
export ANTHROPIC_API_KEY="your_api_key_here"
export OPENAI_API_KEY="your_api_key_here"
```

**On Windows (Command Prompt):**
```cmd
set GEMINI_API_KEY=your_api_key_here
set ANTHROPIC_API_KEY=your_api_key_here
set OPENAI_API_KEY=your_api_key_here
```

**On Windows (PowerShell):**
```powershell
$env:GEMINI_API_KEY="your_api_key_here"
$env:ANTHROPIC_API_KEY="your_api_key_here"
$env:OPENAI_API_KEY="your_api_key_here"
```

### Option 2: .env File

Create a `.env` file in the project root:

```
GEMINI_API_KEY=your_api_key_here
ANTHROPIC_API_KEY=your_api_key_here
OPENAI_API_KEY=your_api_key_here
```

Then load it before running:
```bash
export $(cat .env | xargs)
```

### Getting API Keys

- **Google Gemini** (FREE tier): https://aistudio.google.com/apikey
- **Anthropic Claude**: https://console.anthropic.com/
- **OpenAI GPT**: https://platform.openai.com/api-keys

## Running the Apps

### Variant Extractor (Main App)
```bash
streamlit run app.py --server.port 5000
```

Then open your browser to: `http://localhost:5000`

### Variant Validator
```bash
streamlit run validator_app/app.py --server.port 5000
```

Then open your browser to: `http://localhost:5000`

## Troubleshooting

### Error: "Google Gemini library not available"
**Solution:** Install the package:
```bash
pip install google-genai
```

### Error: "Anthropic library not available"
**Solution:** Install the package:
```bash
pip install anthropic
```

### Error: "OpenAI library not available"
**Solution:** Install the package:
```bash
pip install openai
```

### Error: "GEMINI_API_KEY not set"
**Solution:** Set the environment variable (see "Setting Up API Keys" above)

### Error: "No module named 'pdfplumber'"
**Solution:** Install document parsing libraries:
```bash
pip install pdfplumber openpyxl python-docx
```

### Port Already in Use
If port 5000 is already in use, change to another port:
```bash
streamlit run app.py --server.port 8501
```

## Virtual Environment (Recommended)

To avoid conflicts with other Python projects:

**Create and activate:**
```bash
# Create virtual environment
python -m venv venv

# Activate on macOS/Linux
source venv/bin/activate

# Activate on Windows
venv\Scripts\activate
```

**Install dependencies:**
```bash
pip install streamlit biopython google-genai anthropic openai pandas pdfplumber openpyxl python-docx
```

**Deactivate when done:**
```bash
deactivate
```

## Verification

To verify your installation:

```python
# Run in Python
import streamlit
import biopython
from google import genai
from anthropic import Anthropic
from openai import OpenAI
import pandas
import pdfplumber
import openpyxl
import docx

print("✅ All packages installed successfully!")
```

## Minimum Python Version

Python 3.8 or higher is required.

Check your version:
```bash
python --version
```

## Next Steps

After installation:
1. Set your API key (Gemini recommended for free tier)
2. Run one of the apps
3. Check the README.md for usage instructions

## Support

If you encounter issues:
1. Check this troubleshooting guide
2. Verify Python version (3.8+)
3. Ensure API keys are set correctly
4. Try reinstalling packages in a fresh virtual environment
