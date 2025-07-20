---
tags:
  - Reading
  - Readed
---
# 标题::{{title}}

## 1.Abstract&Info
### 1.1 Abstract
{{abstractNote}}

### 1.2 Info
{% for type, creators in creators | groupby("creatorType") -%}  
{%- for creator in creators -%}  
**{{"First" if loop.first}}{{type | capitalize}}**::
{%- if creator.name %} {{creator.name}}  
{%- else %} {{creator.lastName}}, {{creator.firstName}}  
{%- endif %} 
{% endfor %}~
{%- endfor %}
**Date**:: {{date|format("YYYY")}}
**DOI**: {{DOI}}
**Publication**: {{publicationTitle}}
**PDF**: {{pdfLink}}
**Zotero**: {{pdfZoteroLink}}


## 2. Annotation
{%- macro calloutHeader(type, color) -%}  
{%- if type == "highlight" -%}  
<mark style="background-color: {{color}}">Quote</mark>  
{%- endif -%}  
  
{%- if type == "text" -%}  
Note  
{%- endif -%}  
{%- endmacro -%}

{% persist "annotations" %}
{% set newAnnotations = annotations | filterby("date", "dateafter", lastImportDate) %}
{% if newAnnotations.length > 0 %}
### Imported: {{importDate | format("YYYY-MM-DD h:mm a")}}

{% for each in newAnnotations %}
{{calloutHeader(each.type, each.color)}}
>{{each.annotatedText}} [jump to](zotero://open-pdf/library/items/{{each.attachment.itemKey}}?page={{each.page}}&annotation={{each.id}})

标注：{{each.comment}}
{% endfor%}

{% endif %}
{% endpersist %}

## 3.notes
{% persist “notes” %}


{% endpersist %}

