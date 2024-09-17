# rodent_loss

Tracking the effect of rodents losing 16 miRNA families.

## How to use

You can now create a workspace for each of your analyses using:

```
cd analyses
cookiecutter bionf_cookiecutter --directory="analysis"
```

As soon as an analysis generated a "milestone", 
i.e. figures, tables, or files that will be used in other analyses
save them in the milestones/ directory of the project.

Once you are finished with an analysis, 
add corresponding rules to the snakemake file of the project.

Don't forget to sync this project to GitHub. 
