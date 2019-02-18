# These functions attempt to be colorblind friendly

I am a deuteranope myself (meaning I have a mostly red, but also some green, color vision deficiency.)  So I chose colors that I knew I would be able to tell apart!  I also built in some other methods that make the functions especially colorblindness friendly:

1. The default color palette is built to work for the most common forms of colorblindness.
2. Once the # of colors being used gets too high for the color panel to compensate (8 colors), letters are added by default to try and help.
3. Utilization of shapes, instead of colors/letters, is a built in functionality as well, though this one must be set manually (because combining lettering & shapes does not work well visually, or computationally with ggplot.)
4. I added a **`Simulate()`** function that allows for generation of any my package's plots using a color.panel modified to simulate any of the major types of color blindness.
5. Below, I also describe how to alter your colors for any other package in a way that you can simulate their look to a colorblind individual.

## 4. To simulate looking at any plot generated with this package as if you were colorblind yourself:

Say this is the code you would use to generate the plot:

```
DBDimPlot("age", do.letter=F)
```

The code to visualize this as if you were a deuteranope like me is:

```
Simulate(type = "deutan", plot.function=DBDimPlot, var = "age", do.letter=F)
```

The Simulate() function's inputs are:

- `type` = "deutan", "protan", "tritan" = the type of colorblindness that you want to simulate.  Deutanopia is the most common, and involves primarily red color deficiency, and generally also a bit of green.  Protanopia involves primarily green color deficiency, and generally also a bit of red.  Tritanopia involves primarily blue color deficiency.
- `plot.function` = the function you want to use.  R may try to add (), but delete that if it does.
- `...` = any and all inputs that go into the plotting function you want to use.

## 5. To alter colors going into any plots not generated with my funcitons:

This is based on functions from the colorspace package.  I utilized this package in my Simulate() funciton, as well as the Darken() and Lighten() functions that are also included in the package.  Credit where credit is due and all.

If you have followed the instructions to install this package already, you will have the updated colorspace package.  So:

```
library("colorspace")
deutan.colors <- deutan(colors)
protan.colors <- protan(colors)
tritan.colors <- tritan(colors)
```

then just plug those versions of your colors into whatever function you would like!

#### NOTE: there are varying degrees of colorblindness, but my Simulate() function (and my above code for the deutan(), protan(), and tritan() functions) will simulate for the most severe cases.
