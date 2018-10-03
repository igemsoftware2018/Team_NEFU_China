import qrcode

#text = "39.105.74.74"
text = "http://2018.igem.org/Main_Page"
img = qrcode.make(text)
img.save("./igem.png")
img.show()