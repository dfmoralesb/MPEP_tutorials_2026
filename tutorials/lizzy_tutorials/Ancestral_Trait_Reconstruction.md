# Ancestral State Reconstruction in R
Elizabeth M. Joyce

## Introduction
As discussed during the lecture, there are many different ways to model the evolution of traits. For out practical, we are going to use a simple model to reconstruct the ancestral states of two traits: seed type, and sexual system. Note that these are **categorical traits**, i.e., all of the states are categories (not continuous). There are also multiple options for modelling the evolution of continuous traits, but we will focus on categorical traits for this tutorial.

## Scoring traits

### Seed type in Meliaceae
There are two major seed types in Meliaceae: **winged seeds**, whereby the testa is thin and dry and forms a wing at one end or the whole way around the seed; and seeds without wings (herein referred to as '**unwinged seeds**') whereby the testa is unwinged and either forms a sarcotesta or is covered with an aril.

![](https://github.com/joyceem/MPEP_tutorials/blob/33c6785e67250bf621668255f5ad11c4eb4569af/tutorials/images/Pasted%20image%2020250309105859.png)

A) An example of **winged seeds** in *Toona sinensis* (image: [Roger Culos](https://de.wikipedia.org/wiki/Toona_sinensis#/media/Datei:Toona_sinensis_MHNT.BOT.2010.12.11.jpg)); B) An example of **unwinged seeds** in *Aglaia meridionalis* (image: Elizabeth Joyce), whereby the testa is covered by a bright orange aril.

Did winged seeds evolve once or multiple times? What was the ancestral state of the family? Did winged seeds evolve from unwinged seeds, or vice versa? How frequently do transitions between these states occur? Is it more frequent for winged seeds to evolve into unwinged seeds? Is one state linked to more diversification than the other? These are the sorts of questions that we can tackle with different types of ASR models.

To conduct an ASR, we first need to gather all the information about what we know about seed states in the family. A simple way to do this is by scoring the states of the seed trait for all extant lineages we have in our phylogeny (our tips).

Seeing as there are two potential states, this is a **binary trait**. We will code states of our tips using the following system:
 - 0=unwinged
 - 1= winged
	
Download the `meli_traits_template.xlsx` and use the last generic monograph of the family by Pennington and Styles (1975) to score the states for each of the tips.
### Sexual systems in Meliaceae
There is a great amount of variation in the sexual systems of species in Meliaceae. Members of Meliaceae can be have unisexual flowers in the same individual (monoecious plants); they can have unisexual flowers in distinct individuals (dioecious plants); they can have unisexual and bisexual flowers (polygamous plants), or they can have only bisexual flowers (hermaphroditic plants). This can vary between species within a genus, and sometimes possibly even within a species (but needs further research).

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309113552.png)
Modified from https://commons.wikimedia.org/wiki/File:Monoecy_dioecy_en.svg.

How did these different sexual systems evolve?

Seeing as there are four potential states, this is a **multistate trait**. We will code states of our tips using the following system:
 - 0=unisexual flowers in the same individual (monoecious plants)
 - 1=unisexual flowers in distinct individuals (dioecious plants)
 - 2=unisexual and bisexual flowers (polygamous plants)
 - 3=only bisexual flowers (hermaphroditic plants)
	
Download the `meli_traits_template.xlsx` and use the paper: Laino Gama R, Muellner-Riehl AN, Demarco D, Pirani JR (2021) Evolution of reproductive traits in the mahagony family (Meliaceae). _Journal of Systematics and Evolution_ **59**(1), 21–43. [https://doi.org/10.1111/jse.12572](https://doi.org/10.1111/jse.12572) to score the states for each of the tips.

### Formatting your trait matrix
Different ASR methods will want different formats for your trait matrix. Make sure you understand what the requirements are for the ASR you want to run, and format your table accordingly. I have already done this for you

## Modelling ASR in R
There are different programs you can use to run different types of models (e.g. RevBayes, R, Mesquite, etc). We are going to use R for this tutorial.

Load the instance of [R Studio on the workstation](http://10.153.134.10:8787) by going to your internet browser and typing in http://10.153.134.10:8787 (or clicking the link). Enter your login details.

Set the working directory to  `23_ancestral_trait_reconstruction`. Open the script `ASR_MPEP.R` within R Studio. 

### Evolution of discrete characters with Markov Models
To begin, we will model trait evolution under a simple Markov Model, where the rate class is the same across the whole tree. To do this, we will use the R package `corHMM`. For more information on `corHMM` see the manual for the package: https://cran.r-project.org/web/packages/corHMM/corHMM.pdf, and the paper that describes it: Beaulieu et al. (2013) *Systematic Biology* [doi.org/10.1093/sysbio/syt034](https://doi.org/10.1093/sysbio/syt034).

Read in the chronogram we generated earlier in the week, and the trait data matrix. Check that everything is read in properly, and that all of your taxa in your trait matrix are represented in the tips of your tree.

```
## Read in your tree file

tree <- read.tree("data/Meli_chronogram.nwk")

## Check your tree and the format

tree #tree properties
tree$tip.label 
tree$edge #this is how the tree is stored; relationship between nodes and their descendants
tree$edge.length #this is the length of the branch lengths
plot(tree)

## Now load your trait matrix with states scored for each trait

data <- read.csv("data/meli_traits.csv", header=F) # Matrix of seed type and sexual system

## Check the format of the table

head(data)
  # V1=species
  # V2=seed type: 0=unwinged/1=winged
  # V3=sexual system: 0=unisexual flowers in the same individual (monoecious plants); 1=unisexual flowers in distinct individuals (dioecious plants); 2=unisexual and bisexual flowers (polygamous plants); 3= only bisexual flowers (hermaphroditic plants)
  # polymorphisms handled with "&"
tail(data)
```
#### Modelling seed evolution (binary trait)
Now, we are going to compare the results of an ER (Equal Rates) vs. ARD (All Rates Different) model of evolution for our binary seed trait.

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309120446.png)

Run your ER model:
```
ans <- rayDISC(tree, data, ntraits=1, charnum=1, model="ER")
```

Look at the results of the model. 
```
> ans
Fit
      -lnL      AIC     AICc ntax
 -11.05538 24.11076 24.27076   27

Rates
            0           1
0          NA 0.004812369
1 0.004812369          NA

Arrived at a reliable solution
```
Note the equal rates in the transition rate matrix: this is what we asked for with the ER model!

Once you've followed the script to check out the format of the ancestral states, plot your results. You should see something like this:

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309121027.png)

Now, if we want to get the values for the proportional likelihood (~probability) of each ancestral state, we need to identify the node number, and extract the values for that node.

```
melioideae <- getMRCA(tree, tip=c("MELI_Quivisianthe_papinae", "MELI_Pterorhachis_zenkeri"))

melioideae #This is a node/edge number in ape's language, to get the corresponding row in rayDISC's ans$states, we need to substract the number of tips

nrtips <- length(tree$tip.label)
ans$states[melioideae-nrtips,]
```
What is the most likely seed type in the ancestor Melioideae according to this model? 

Can you do the same for the ancesral node of the Cedreloideae subfamily?

Now, we are going to model evolution again according to a Markov Model, but this time we will use an ARD transition rate matrix.

```
ans <- rayDISC(tree, data, ntraits=1, charnum=1, model="ARD")
```

Note that this time, we have unequal rates in the transition rate matrix (this is what we asked for with the ARD model):
```
> ans
Fit
      -lnL      AIC     AICc ntax
 -10.36366 24.72732 25.22732   27

Rates
            0           1
0          NA 0.002967107
1 0.008803931          NA

Arrived at a reliable solution
```
What does this suggest about how seeds evolve in Meliaceae?

Plot the results of the ARD model.

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309122019.png)

How do the results differ from the ER model?

But which one is a better model for our data? We can test this by comparing the AIC:
```
> ansER$AIC
[1] 24.11076
> ansARD$AIC
[1] 24.72732
> deltaAIC <- ansER$AIC - ansARD$AIC
> deltaAIC
[1] -0.6165636
```
We can see that the ER is a slightly better fit than the ARD model.

#### Modelling sexual system evolution (multistate trait)
Follow the R script to repeat what we just did for seeds to reconstruct the ancestral state of sexual systems in Meliaceae.

When you plot your results, you should see something like this for the ER...

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309124021.png)

...and this for the ARD:

![](https://github.com/joyceem/MPEP_tutorials/blob/90a90e886ffa5f455e03c8b6015a51c2c8228c94/tutorials/images/Pasted%20image%2020250309124131.png)

Explore the output. Which is the best model? What are the most likely ancestral states for Melioideae, Cedreloideae and Meliaceae?

### Evolution of discrete characters with rate changes
One potential problem with traditional Mk models such as the ones we just ran is that they assume that transition rates are fixed across the whole tree. One way we can allow for rates to change across the tree is by modelling states that allow for different rate categories across the tree. One optoin for this is by applying what's called a "Hidden Markov Model", whereby we add an extra "hidden" state that represent different rate categories. Note that this is different from Hidden SSE models, which we will move on to later, as we are not modelling anything about speciation or extinction rates. Let's go back to our seed data to try this.

First, let's rerun an ER model with one rate category:
```
ansER_HMM1 <-corHMM(phy = tree, data = seeds_data, rate.cat = 1, model="ER")
ansER_HMM1
```
Inspect the output. You should see something very similar to out previous ER model, as it is the same model:
```
> ansER_HMM1
Fit
      -lnL      AIC     AICc Rate.cat ntax
 -11.05538 24.11076 24.27076        1   27

Legend
  1   2 
"0" "1" 

Rates
            (1,R1)      (2,R1)
(1,R1)          NA 0.004812382
(2,R1) 0.004812382          NA

Arrived at a reliable solution
```

Now, allow for two rate categories (fast=R1 and slow=R2).
```
ansER_HMM2 <-corHMM(phy = tree, data = seeds_data, rate.cat = 2, model="ER")
ansER_HMM2
```
Inspect the output. We can see that even though we have allowed for two rate categories, most transitions are still happening in the same rate class:
```
> ansER_HMM2$solution
             (1,R1)       (2,R1)       (1,R2)       (2,R2)
(1,R1)           NA  0.004812301 2.194261e-09           NA
(2,R1)  0.004812301           NA           NA 2.194261e-09
(1,R2) 78.431733784           NA           NA 4.249587e-06
(2,R2)           NA 78.431733784 4.249587e-06           NA
```

Now repeat for three rate categories.

Compare the AIC of all models. Which one is the most likely?
```
> ansER_HMM1$AIC
[1] 24.11076
> ansER_HMM2$AIC
[1] 30.11076
> ansER_HMM3$AIC
[1] 39.31278
```

We can see that even though we have allowed for the rate categories to differ across the tree, the model with only one category is still the best. This suggests that the rate of seed evolution is constant across the whole tree; but this might also be due to having such a small phylogeny. The results could be different if we had more tips. What results might we expect to see if Cedreloideae had a faster rate of seed evolution than Meiloideae?

