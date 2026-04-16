library(tidyverse)
library(GGally)
library(network)
library(sna)


# random graph
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:10]

# plot
net %v% "phono" = ifelse(letters[1:10] %in% c("b", "e", "i", "m", "x"), "vowel", "consonant")
net %v% "color" = ifelse(net %v% "phono" == "vowel", "steelblue", "tomato")
p <- ggnet2(net, color = "color", size = 2, ) +
  theme_void() + theme_transparent()


# PNG hex sticker
sticker(p, 
        package="PolyLinkR", p_size = 7.5, p_color = "tomato", p_y = 0.55,
        s_x = 1, s_y = 1.2,
        s_width = 1, s_height = 1,
        h_fill = "grey15", h_color = "steelblue",
        filename= file.path("inst/sticker/polylinkr.png"))

# SVG hex sticker
sticker(p, 
        package="PolyLinkR", p_size = 7.5, p_color = "tomato", p_y = 0.55,
        s_x = 1, s_y = 1.2,
        s_width = 1, s_height = 1,
        h_fill = "grey15", h_color = "steelblue",
        filename= file.path("inst/sticker/polylinkr.svg"))


