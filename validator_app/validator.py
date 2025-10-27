import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pubmed_api import PubMedAPI
import json

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
    from google.genai import types
    HAS_GEMINI = True
except ImportError:
    HAS_GEMINI = False

class VariantValidator:
    def __init__(self, api_key: str, model_choice: str):
        self.api_key = api_key
        self.model_choice = model_choice
        self.pubmed = PubMedAPI()
        
        if "Gemini" in model_choice:
            if not HAS_GEMINI:
                raise Exception("Google Gemini library not available")
            self.client = genai.Client(api_key=api_key)
            self.model = "gemini-2.5-flash"
            self.llm_type = "gemini"
        elif "Claude" in model_choice:
            if not HAS_ANTHROPIC:
                raise Exception("Anthropic library not available")
            self.client = Anthropic(api_key=api_key)
            self.model = "claude-3-haiku-20240307"
            self.llm_type = "anthropic"
        else:
            if not HAS_OPENAI:
                raise Exception("OpenAI library not available")
            self.client = OpenAI(api_key=api_key)
            self.model = "gpt-4o-mini"
            self.llm_type = "openai"
    
    def call_llm(self, prompt: str) -> str:
        """Call the selected LLM with a prompt"""
        try:
            if self.llm_type == "anthropic":
                message = self.client.messages.create(
                    model=self.model,
                    max_tokens=2048,
                    messages=[{"role": "user", "content": prompt}]
                )
                return message.content[0].text
            elif self.llm_type == "openai":
                response = self.client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": "You are a biomedical research assistant. Always respond with valid JSON only."},
                        {"role": "user", "content": prompt}
                    ],
                    response_format={"type": "json_object"},
                    max_tokens=2048
                )
                return response.choices[0].message.content
            else:
                response = self.client.models.generate_content(
                    model=self.model,
                    contents=[types.Content(role="user", parts=[types.Part(text=prompt)])],
                    config=types.GenerateContentConfig(
                        response_mime_type="application/json"
                    )
                )
                return response.text
        except Exception as e:
            raise Exception(f"LLM call failed: {str(e)}")
    
    def validate_variant(self, pmid: str, annotated_variant: str, annotated_phenotype: str, gene: str = "") -> dict:
        """
        Validate a variant annotation against the source paper
        
        Returns:
            dict with keys: validation_status, validation_reason, paper_title, reported_data
        """
        try:
            pmid_list = [pmid]
            
            import re
            gene_from_variant = gene if gene else "ValidationQuery"
            
            if not gene and annotated_variant:
                variant_parts = annotated_variant.split('|')
                for part in variant_parts:
                    part = part.strip()
                    gene_match = re.match(r'^([A-Z][A-Z0-9]{2,})\b', part)
                    if gene_match:
                        gene_from_variant = gene_match.group(1)
                        break
            
            articles = self.pubmed.fetch_article_details(pmid_list, gene=gene_from_variant, variant=None)
            
            if not articles:
                return {
                    'pmid': pmid,
                    'validation_status': 'UNSURE',
                    'validation_reason': 'Paper not found in PubMed',
                    'paper_title': '',
                    'reported_data': ''
                }
            
            article = articles[0]
            paper_title = article.get('title', '')
            
            text_content = f"Title: {paper_title}\n\nAbstract: {article.get('abstract', '')}"
            
            full_text = article.get('full_text', '')
            if full_text:
                text_content += f"\n\nFull Text:\n{full_text[:20000]}"
            
            validation_prompt = f"""You are validating genetic variant annotations against the source publication.

**Annotated Information from CSV:**
- Variant: {annotated_variant}
- Phenotype/Diagnosis: {annotated_phenotype}

**Source Paper Content:**
{text_content[:50000]}

**Task:**
Compare the annotated information against what is actually reported in the paper. Determine:
1. Is the variant correctly identified? (check variant name, HGVS notation, gene)
2. Is the phenotype/diagnosis correctly annotated? (check reported clinical features, disease, outcomes)
3. Overall assessment: CORRECT, INCORRECT, or UNSURE

**Output JSON format:**
{{
    "validation_status": "CORRECT" | "INCORRECT" | "UNSURE",
    "variant_match": "exact" | "partial" | "mismatch" | "not_found",
    "phenotype_match": "exact" | "partial" | "mismatch" | "not_found",
    "reported_variant": "what the paper actually says about the variant",
    "reported_phenotype": "what the paper actually says about phenotype/diagnosis",
    "validation_reason": "detailed explanation of your assessment"
}}

**Validation Guidelines:**
- CORRECT: Both variant and phenotype match the paper (exact or semantically equivalent)
- INCORRECT: Clear mismatch between annotation and paper (wrong variant, wrong phenotype, or contradictory info)
- UNSURE: Paper doesn't provide enough information, ambiguous reporting, or unclear match
- Be strict but fair - minor notation differences (p. vs protein name) are OK if semantically identical
- Consider that abbreviations and full terms may be used interchangeably
"""

            response = self.call_llm(validation_prompt)
            
            try:
                result_data = json.loads(response)
            except json.JSONDecodeError:
                json_match = response.find('{')
                if json_match != -1:
                    json_end = response.rfind('}') + 1
                    result_data = json.loads(response[json_match:json_end])
                else:
                    result_data = {
                        'validation_status': 'UNSURE',
                        'validation_reason': 'AI response parsing error',
                        'reported_variant': '',
                        'reported_phenotype': ''
                    }
            
            return {
                'pmid': pmid,
                'validation_status': result_data.get('validation_status', 'UNSURE'),
                'validation_reason': result_data.get('validation_reason', 'No reason provided'),
                'paper_title': paper_title,
                'reported_data': f"Variant: {result_data.get('reported_variant', 'N/A')} | Phenotype: {result_data.get('reported_phenotype', 'N/A')}"
            }
            
        except Exception as e:
            return {
                'pmid': pmid,
                'validation_status': 'ERROR',
                'validation_reason': f'Validation error: {str(e)}',
                'paper_title': '',
                'reported_data': ''
            }
