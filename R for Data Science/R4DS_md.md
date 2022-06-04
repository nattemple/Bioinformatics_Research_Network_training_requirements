---
title: "R4DS Review, Nathan Temple"
output: html_notebook
editor_options: 
  chunk_output_type: inline
---

**The R script**
Thank you for quickly solving the issues. The report is pretty clear now. However, there are still two major problems.

1. **Population density**, you found Singapore has the higher population density. However, the correct answer is not Singapore. I searched for the correct answer within the table you provided. However, the correct answer is not there. You seem to somehow filter that value.

I had run 'drop_na' on the dataset while including continent, which filtered out Macao SAR, China and Monaco.  After removing continent and re-running, these 2 countries are tied for the highest average rank.

2. The same issue goes for the **life expectancy** question. I suggest you use the raw **gapminder** data rather than the processed version for these questions.

I had the same issue of running 'drop_na' while including continent.  After dropping continent and re-running, Maldives had the highest change in life expectancy over the period.

**Style**

1. I highly recommend using  *code_folding: hide* so we can see your code chunks if we want.

https://bookdown.org/yihui/rmarkdown-cookbook/fold-show.html

```
 ---
title: "R4DS"
author: "Nathan Temple"
Date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: hide
--- 
  
```
I've included the code_floating option in my Rmd.  I'm not sharing any code in this assignment but will use it in the next training project (RNA-Seq Analysis).  I also included a table of contents.


2. If you want to plot two figures for the same question you can use **tabs**. 

https://bookdown.org/yihui/rmarkdown-cookbook/html-tabs.html

I switched to using tabs for three pairs of charts.

3. Please try to avoid printing unfriendly elements while printing **datatables**.For example, the first column "x1" seems to be a mistake and the size of the table is larger than the screen.
https://stackoverflow.com/questions/30765338/how-to-make-the-horizontal-scrollbar-visible-in-dtdatatable

I removed the datatable since it's too much information for the purposes of this project.  I replaced it with a brief narrative description of the source data.

4. The axis titles on your plotly plot are also unfriendly and should be changed.

I modified the axis titles to be more readable.

Keep up the good work! **- BRN SA Team**

