import os
import json
import sys
from typing import Dict, List, Optional

try:
    from anthropic import Anthropic
    HAS_ANTHROPIC = True
except ImportError:
    HAS_ANTHROPIC = False

try:
    from openai import OpenAI
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False

try:
    import google.generativeai as genai
    from google.generativeai.types import GenerationConfig
    HAS_GEMINI = True
except ImportError:
    HAS_GEMINI = False

class LLMAnalyzer:
    def __init__(self, model_choice: str = "anthropic"):
        self.model_choice = model_choice
        
        if model_choice == "anthropic":
            if not HAS_ANTHROPIC:
                raise Exception("Anthropic library not available")
            anthropic_key = os.environ.get('ANTHROPIC_API_KEY')
            if not anthropic_key:
                raise Exception("ANTHROPIC_API_KEY not set")
            self.client = Anthropic(api_key=anthropic_key)
            self.model = "claude-3-haiku-20240307"
        elif model_choice == "openai":
            if not HAS_OPENAI:
                raise Exception("OpenAI library not available")
            openai_key = os.environ.get('OPENAI_API_KEY')
            if not openai_key:
                raise Exception("OPENAI_API_KEY not set")
            self.client = OpenAI(api_key=openai_key)
            self.model = "gpt-3.5-turbo"
        elif model_choice == "gemini":
            if not HAS_GEMINI:
                raise Exception("Google Gemini library not available")
            gemini_key = os.environ.get('GEMINI_API_KEY')
            if not gemini_key:
                raise Exception("GEMINI_API_KEY not set")
            genai.configure(api_key=gemini_key)
            self.model = genai.GenerativeModel("gemini-1.5-flash")
        else:
            raise Exception(f"Unsupported model choice: {model_choice}")
    
    def extract_variant_data(
        self,
        article_text: str,
        gene: str,
        variant: Optional[str] = None,
        target_phenotypes: Optional[List[str]] = None
    ) -> Dict:
        phenotype_focus = "None specified"
        if target_phenotypes:
            phenotype_focus = ", ".join(target_phenotypes)
        prompt = f"""Analyze the following research article and extract genetic variant and clinical data at TWO levels:

Gene of interest: {gene}
{f'Variant of interest: {variant}' if variant else ''}

Prioritize extracting or confirming the following phenotypes when available: {phenotype_focus}.

**LEVEL 1: INDIVIDUAL PATIENT DATA** (case reports, patient narratives)
For each individual patient/participant mentioned, extract:
1. Genetic variant (gene name and specific variant, like p.Tyr54Asn or Y54N)
2. Genotype (heterozygous, homozygous, compound heterozygous)
3. Phenotypes/diagnoses
4. Age, sex, ethnicity/ancestry (if mentioned)
5. Treatment and clinical outcomes (if mentioned)

**LEVEL 2: AGGREGATE STUDY DATA** (cohort studies, population genetics, meta-analyses)
Extract study-level statistics from tables and supplemental materials:
1. Study design (case-control, cohort, GWAS, meta-analysis, family study, etc.)
2. Sample sizes (number of cases with variant, number of controls, total cohort size)
3. Ethnicity/ancestry breakdown (e.g., "500 European, 200 East Asian, 100 African")
4. Study location/geography (country, institution, population source)
5. Cohort name if mentioned (e.g., "UK Biobank", "CHARGE Consortium", "Framingham Heart Study")
6. Variant frequency in cases vs controls (if provided)
7. For tables: extract ALL variants in the gene of interest with their case/control counts

**DEDUPLICATION METADATA**
Extract information to help identify potential participant overlap:
- Cohort/study names
- Recruitment dates/periods
- Geographic location + institution
- Author names/institutions (for tracking follow-up studies)

Article text:
{article_text[:50000]}

Respond ONLY with valid JSON in this exact format:
{{
  "has_variant_data": true/false,
  "study_metadata": {{
    "study_type": "case-control/cohort/GWAS/meta-analysis/case-report/family-study/unknown",
    "location": "country, city, or institution if mentioned, else null",
    "cohort_name": "official cohort/study name if mentioned, else null",
    "recruitment_period": "dates if mentioned, else null",
    "primary_institution": "lead institution/affiliation if mentioned, else null"
  }},
  "individual_patients": [
    {{
      "gene": "gene name",
      "variant": "variant notation",
      "genotype": "heterozygous/homozygous/compound het/unknown",
      "phenotypes": ["phenotype1", "phenotype2"],
      "age": "age if mentioned or null",
      "sex": "sex if mentioned or null",
      "ethnicity": "ethnicity/ancestry if mentioned or null",
      "treatment": "treatment if mentioned or null",
      "outcome": "outcome if mentioned or null"
    }}
  ],
  "aggregate_data": [
    {{
      "gene": "gene name",
      "variant": "variant notation",
      "cases_with_variant": "number or null",
      "total_cases": "number or null",
      "controls_with_variant": "number or null",
      "total_controls": "number or null",
      "ethnicity_breakdown": "description of population ancestry, e.g. '500 EUR, 200 EAS' or null",
      "allele_frequency_cases": "frequency if mentioned or null",
      "allele_frequency_controls": "frequency if mentioned or null",
      "phenotype": "associated phenotype/condition",
      "notes": "any relevant context from table or text"
    }}
  ],
  "confidence": "high/medium/low"
}}

If no variant data is found, return: {{"has_variant_data": false, "study_metadata": {{}}, "individual_patients": [], "aggregate_data": [], "confidence": "low"}}
"""

        try:
            if self.model_choice == "anthropic":
                message = self.client.messages.create(
                    model=self.model,
                    max_tokens=2048,
                    messages=[
                        {"role": "user", "content": prompt}
                    ]
                )
                response_text = message.content[0].text
            elif self.model_choice == "openai":
                response = self.client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": "You are a biomedical research assistant. Always respond with valid JSON only."},
                        {"role": "user", "content": prompt}
                    ],
                    response_format={"type": "json_object"},
                    max_tokens=2048
                )
                response_text = response.choices[0].message.content
            else:
                response = self.model.generate_content(
                    prompt,
                    generation_config=GenerationConfig(
                        response_mime_type="application/json"
                    )
                )
                response_text = response.text
            
            cleaned_response = self._clean_json_response(response_text)
            result = json.loads(cleaned_response)
            return result

        except json.JSONDecodeError as e:
            return {
                "has_variant_data": False,
                "study_metadata": {},
                "individual_patients": [],
                "aggregate_data": [],
                "confidence": "low",
                "error": f"JSON parsing error: {str(e)}"
            }
        except Exception as e:
            return {
                "has_variant_data": False,
                "study_metadata": {},
                "individual_patients": [],
                "aggregate_data": [],
                "confidence": "low",
                "error": str(e)
            }

    def batch_analyze_articles(
        self,
        articles: List[Dict],
        gene: str,
        variant: Optional[str] = None,
        target_phenotypes: Optional[List[str]] = None
    ) -> List[Dict]:
        results = []

        for article in articles:
            analysis = self.extract_variant_data(
                article['full_text'],
                gene,
                variant,
                target_phenotypes=target_phenotypes
            )

            result = {
                "pmid": article['pmid'],
                "title": article['title'],
                "year": article['year'],
                "authors": article['authors'],
                "has_variant_data": analysis.get('has_variant_data', False),
                "study_metadata": analysis.get('study_metadata', {}),
                "individual_patients": analysis.get('individual_patients', []),
                "aggregate_data": analysis.get('aggregate_data', []),
                "confidence": analysis.get('confidence', 'low'),
                "error": analysis.get('error', None)
            }

            results.append(result)

        return results

    @staticmethod
    def _clean_json_response(response_text: str) -> str:
        """Attempt to extract valid JSON from an LLM response."""
        text = response_text.strip()

        if text.startswith("```"):
            parts = [segment.strip() for segment in text.split("```") if segment.strip()]
            for segment in parts:
                if "{" in segment:
                    text = segment
                    break
            else:
                text = parts[-1] if parts else text

            # remove language hints like ```json
            if text.lower().startswith("json"):
                text = text[4:].strip()

        first_brace = text.find("{")
        last_brace = text.rfind("}")
        if first_brace != -1 and last_brace != -1:
            text = text[first_brace:last_brace + 1]

        return text
