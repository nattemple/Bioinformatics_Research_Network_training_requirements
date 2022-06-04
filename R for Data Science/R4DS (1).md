---
title: "R4DS Review, Nathan Temple"
output: html_notebook
---

**The R script**

1. You used pipes *%>%* for some of the operations while some of them you didn't. Why not use pipe for all?

I switched several assignment (<-) statements to pipe (%>%).  In some cases, I retained assignment statements (<-) to retain intermediate datasets for downstream processing or browsing while developing.

2. The *print* operation is not necessary. You can directly use **plot + ggtitle...**
  ```
  print(Plot + ggtitle("CO2 emissions vs. GDP") +
  labs(x = "CO2 emission (metric tons) per Cap", y = "GDP Per Cap") +
  theme(plot.title = element_text(hjust = 0.5)))
  ```
I am no longer using *print* statements.

3. Because of the one outlier point, the whole plot seems to be skewed. What can you do to have a better looking and,interpretable plot? **Hint:** What are your reflections on this plot? Do you need data transformation? Or, can you set set boundaries for each axis? 

I added a reflection that there is a strong positive correlation between the variables.  I also added a second plot with the outlier removed so that readers can review the remaining points in more detail.

4. *Correlation between`'CO2 emissions (metric tons per capita)'` and `gdpPercap`* You can round the values of the correlations. I.e. instead of **0.9387918**, **0.939** can be more interpretable. Additionally, you choose to apply *Pearson correlation*. Can you justify your choice? What are the statistical reasons behind your choice?

I rounded the correlations to 3 decimal digits.  I also plotted 'CO2 emissions (metric tons per capita)'` and `gdpPercap` individually and found that they do not have a normally shaped distribution, so I switched correlations to Spearman.  The 3 assumptions of the Spearman correlation are satisfied, since the variables are interval, the variables are paired (they are 2 quantities drawn from a country in a specific year), and there is a monotonic relationship between the variables (as one variable increases, the second variable generally increases).
 
5. *CO2 Emissions vs. GDP interactive Plot* Again, your plot seems to be pressed. The point of using *plotly* is to interactively interpret data points. However, since all the data has collided on one side, that's almost impossible.

Switching to the Spearman correlation changed the year I filtered to (from 1967 to 2002) and resulted in a better spread in the dataset and plot.

6. *Relationship between Continent and Energy Use* As mentioned, a statistical test is required to test the relationship between *continent* and *Energy use (kg of oil equivalent per capita)*. You applied **one-way ANOVA test** for the data filtered for the year 1987. Why not apply it to all data? Then, you can plot your findings through box plots via faceting through the years. 

I've plotted the data for all years and also added a faceted plot for 3 of the years as a visual sample.  The Anova test previously included all years and I've retained that in the current version.

8. *Imports of Goods and Service in Europe vs. Asia* `An ANOVA test, signifcant only at a p-value of 0.158, does not allow us to reject the null hypothesis that the Imports of Europe and Asia during this time period are significantly different:`
What is your null hypothesis? What is it that you reject? Did you use **one-way ANOVA** or **two-way ANOVA**? In addition, you don't need the first plot, the box plot you have is enough to justify your findings.

I've added more detail on the one-way Anova and the null hypothesis to the narrative.  I retained the scatter plot as it provides some usefuly contextual information on the data.

9. *Population Density*  Did you calculate the average population density per country? Your answer to this question is incorrect.

I clarified the meaning of this question with Henry.  I calculated the within-year rank (of population density) for each country and then averaged the rank for each country across all years.  I found that Singapore still had the highest average ranking across all years.

10. *Life Expectancy* What is delta? What does it mean, higher delta? The printed table seems meaningless. Can you make a plot? Which country has shown the greatest increase in life expectancy since 1962? *Formula*: life expectancy = the value from 2007 minus the value from 1962. You should calculate the expectancy then, plot or make a table from the results in order.

I re-labelled the columns and added commentary that the table is ordered from highest to lowest life-expectancy change for more clarity.


**Style**

1. You should load the libraries at the beginning, not in the middle.
  ```
library(plotly)

I moved all library statements to the beginning.
  ```
2. Should've explained what you did for which task. Instead of using one code chunk, you can use chunks and labels of which tasks that you're trying to answer.

I labelled each section with the question number for more clarity.  I've also used dataset names to indicate the purpose of each data processing step.

3. Instead of printing tables, you can use *datatable* package to have more good looking tables.

I've switched to using *datatable* for displaying table in the Rmd.


**Final Notes:** Overall, the goal of this training is to introduce exploratory data analysis or, in the way 
I call it, data-based storytelling. We are focusing on how you represent and, express your findings. Even though it is not required to have a perfectly styled **markdown** for this training, I highly encourage you to do so. Therefore, I am attaching these sources for you.

https://google.github.io/styleguide/docguide/style.html
https://about.gitlab.com/handbook/markdown-guide/
https://bookdown.org/yihui/rmarkdown-cookbook/chunk-styling.html
tidyverse.org/blog/2017/12/styler-1.0.0/

I used Lintr and Styler on the code.

Good luck and, keep up the good work! **- BRN SA Team**
