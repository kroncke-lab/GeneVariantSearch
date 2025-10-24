import os
import time
from Bio import Entrez
import xml.etree.ElementTree as ET
from typing import List, Dict, Optional

Entrez.email = os.environ.get("PUBMED_EMAIL", "user@example.com")

class PubMedAPI:
    def __init__(self, email: str = None):
        if email:
            Entrez.email = email
        self.max_results = 20
    
    def search_publications(self, gene: str, variant: Optional[str] = None, max_results: int = 20) -> List[str]:
        query = gene
        if variant:
            query = f"{gene} AND {variant}"
        
        query += " AND (genetic variant OR mutation OR polymorphism OR genotype)"
        
        try:
            search_handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort="relevance"
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            id_list = search_results.get("IdList", [])
            return id_list
        except Exception as e:
            raise Exception(f"PubMed search failed: {str(e)}")
    
    def get_pmc_id(self, pmid: str) -> Optional[str]:
        """Get PMC ID for a given PMID if available in PubMed Central"""
        try:
            link_handle = Entrez.elink(
                dbfrom="pubmed",
                db="pmc",
                id=pmid,
                linkname="pubmed_pmc"
            )
            link_results = Entrez.read(link_handle)
            link_handle.close()
            
            if link_results and link_results[0].get("LinkSetDb"):
                pmc_ids = link_results[0]["LinkSetDb"][0]["Link"]
                if pmc_ids:
                    return pmc_ids[0]["Id"]
            return None
        except:
            return None
    
    def fetch_pmc_fulltext(self, pmc_id: str) -> Optional[str]:
        """Fetch full text from PMC including supplemental data"""
        try:
            fetch_handle = Entrez.efetch(
                db="pmc",
                id=pmc_id,
                rettype="xml",
                retmode="xml"
            )
            xml_data = fetch_handle.read()
            fetch_handle.close()
            
            root = ET.fromstring(xml_data)
            
            full_text_parts = []
            
            title = root.find(".//article-title")
            if title is not None and title.text:
                full_text_parts.append(f"TITLE: {title.text}")
            
            abstract = root.find(".//abstract")
            if abstract is not None:
                abstract_text = " ".join([elem.text or "" for elem in abstract.iter() if elem.text])
                if abstract_text.strip():
                    full_text_parts.append(f"ABSTRACT: {abstract_text}")
            
            body = root.find(".//body")
            if body is not None:
                body_text = " ".join([elem.text or "" for elem in body.iter() if elem.text and elem.tag in ['p', 'title', 'td', 'th']])
                if body_text.strip():
                    full_text_parts.append(f"MAIN TEXT: {body_text[:15000]}")
            
            tables = root.findall(".//table-wrap")
            if tables:
                table_texts = []
                for idx, table in enumerate(tables[:10]):
                    table_content = " ".join([elem.text or "" for elem in table.iter() if elem.text])
                    if table_content.strip():
                        table_texts.append(f"Table {idx+1}: {table_content}")
                if table_texts:
                    full_text_parts.append(f"TABLES: {' | '.join(table_texts)}")
            
            supplementary = root.findall(".//supplementary-material")
            if supplementary:
                supp_texts = []
                for supp in supplementary:
                    supp_content = " ".join([elem.text or "" for elem in supp.iter() if elem.text])
                    if supp_content.strip():
                        supp_texts.append(supp_content)
                if supp_texts:
                    full_text_parts.append(f"SUPPLEMENTARY: {' | '.join(supp_texts)}")
            
            if full_text_parts:
                return "\n\n".join(full_text_parts)
            return None
            
        except Exception as e:
            return None
    
    def fetch_article_details(self, pmid_list: List[str]) -> List[Dict]:
        if not pmid_list:
            return []
        
        articles = []
        
        try:
            fetch_handle = Entrez.efetch(
                db="pubmed",
                id=pmid_list,
                rettype="xml",
                retmode="xml"
            )
            
            xml_data = fetch_handle.read()
            fetch_handle.close()
            
            root = ET.fromstring(xml_data)
            
            for article in root.findall(".//PubmedArticle"):
                try:
                    pmid_elem = article.find(".//PMID")
                    pmid = pmid_elem.text if pmid_elem is not None else "Unknown"
                    
                    title_elem = article.find(".//ArticleTitle")
                    title = title_elem.text if title_elem is not None else "No title"
                    
                    abstract_texts = article.findall(".//AbstractText")
                    abstract = " ".join([
                        elem.text for elem in abstract_texts if elem.text
                    ]) if abstract_texts else "No abstract available"
                    
                    year_elem = article.find(".//PubDate/Year")
                    year = year_elem.text if year_elem is not None else "Unknown"
                    
                    authors = []
                    for author in article.findall(".//Author"):
                        lastname = author.find("LastName")
                        forename = author.find("ForeName")
                        if lastname is not None and forename is not None:
                            authors.append(f"{forename.text} {lastname.text}")
                    
                    author_str = ", ".join(authors[:3]) if authors else "Unknown authors"
                    if len(authors) > 3:
                        author_str += " et al."
                    
                    pmc_id = self.get_pmc_id(pmid)
                    content_type = "abstract"
                    full_text = f"{title}\n\n{abstract}"
                    
                    if pmc_id:
                        time.sleep(0.34)
                        pmc_fulltext = self.fetch_pmc_fulltext(pmc_id)
                        if pmc_fulltext:
                            full_text = pmc_fulltext
                            content_type = "full_text_with_supplements"
                    
                    articles.append({
                        "pmid": pmid,
                        "title": title,
                        "abstract": abstract,
                        "year": year,
                        "authors": author_str,
                        "full_text": full_text,
                        "content_type": content_type,
                        "pmc_id": pmc_id
                    })
                    
                except Exception as e:
                    continue
            
            return articles
            
        except Exception as e:
            raise Exception(f"Failed to fetch article details: {str(e)}")
    
    def search_and_fetch(self, gene: str, variant: Optional[str] = None, max_results: int = 20) -> List[Dict]:
        pmid_list = self.search_publications(gene, variant, max_results)
        
        if not pmid_list:
            return []
        
        time.sleep(0.34)
        
        articles = self.fetch_article_details(pmid_list)
        return articles
