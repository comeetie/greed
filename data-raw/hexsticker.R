library(hexSticker)
imgurl <- "https://www.comeetie.fr/greed.jpg"
library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Fredoka One", "fredoka")
## Automatically use showtext to render text for future devices
showtext_auto()
sticker(imgurl, package="GREED", p_size=22, s_x=1, s_y=.95, s_width=.6,p_y=1.5,
        filename="greed.png",h_color="#fe0096",h_fill="#ffffff",p_color="#000000",p_family = "fredoka")
