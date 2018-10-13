import qrcode

#text = "39.105.74.74"
text = "http://2018.igem.org/Team:NEFU_China"
img = qrcode.make(text)
img.save("./igem_nefu.png")
img.show()