##### PACKAGES #####

library(tidyverse)
library(readxl)
library(ggpubr)
library(scales)
library(extrafont)
library(ggrepel)
library(ggh4x)
library(igraph)
library(ggraph)
library(ggtext)
library(tidygraph)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(base::union)

sample(x = (1:1000), size = 1000, replace = T)

