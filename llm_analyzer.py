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
    from google import genai
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
            self.client = genai.Client(api_key=gemini_key)
            self.model = "gemini-2.5-flash"
        else:
            raise Exception(f"Unsupported model choice: {model_choice}")
    
    def extract_variant_data(self, article_text: str, gene: str, variant: str = None) -> Dict:
        prompt = f"""Analyze the following research article and extract genetic variant and clinical data.

Gene of interest: {gene}
{f'Variant of interest: {variant}' if variant else ''}

Extract the following information if present:
1. Genetic variants mentioned (gene name and specific variant notation like p.Tyr54Asn)
2. Carrier genotypes (heterozygous, homozygous, compound heterozygous)
3. Clinical phenotypes (symptoms, conditions, measurements like QT prolongation, arrhythmia, syncope)
4. Patient demographics (age, sex) if mentioned
5. Treatment information if mentioned
6. Clinical outcomes if mentioned

Article text:
{article_text[:8000]}

Respond ONLY with valid JSON in this exact format:
{{
  "has_variant_data": true/false,
  "variants": [
    {{
      "gene": "gene name",
      "variant": "variant notation",
      "genotype": "heterozygous/homozygous/compound het/unknown",
      "phenotypes": ["phenotype1", "phenotype2"],
      "age": "age if mentioned or null",
      "sex": "sex if mentioned or null",
      "treatment": "treatment if mentioned or null",
      "outcome": "outcome if mentioned or null"
    }}
  ],
  "confidence": "high/medium/low"
}}

If no variant data is found, return: {{"has_variant_data": false, "variants": [], "confidence": "low"}}
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
                response = self.client.models.generate_content(
                    model=self.model,
                    contents=prompt
                )
                response_text = response.text
            
            result = json.loads(response_text)
            return result
            
        except json.JSONDecodeError as e:
            return {
                "has_variant_data": False,
                "variants": [],
                "confidence": "low",
                "error": f"JSON parsing error: {str(e)}"
            }
        except Exception as e:
            return {
                "has_variant_data": False,
                "variants": [],
                "confidence": "low",
                "error": str(e)
            }
    
    def batch_analyze_articles(self, articles: List[Dict], gene: str, variant: str = None) -> List[Dict]:
        results = []
        
        for article in articles:
            analysis = self.extract_variant_data(article['full_text'], gene, variant)
            
            result = {
                "pmid": article['pmid'],
                "title": article['title'],
                "year": article['year'],
                "authors": article['authors'],
                "has_variant_data": analysis.get('has_variant_data', False),
                "variants": analysis.get('variants', []),
                "confidence": analysis.get('confidence', 'low'),
                "error": analysis.get('error', None)
            }
            
            results.append(result)
        
        return results
