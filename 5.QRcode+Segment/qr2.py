from PIL import Image
import qrcode, os


def create_qrcode(url, qrcodename):
    qr = qrcode.QRCode(
        version=1,
        error_correction=qrcode.ERROR_CORRECT_H,
        box_size=8,
        border=1,
    )

    qr.add_data(url)
    qr.make(fit=True)

    img = qr.make_image()
    img = img.convert("RGBA")
    logo = Image.open('logo.png')

    w, h = img.size
    logo_w, logo_h = logo.size
    factor = 4
    s_w = int(w / factor)
    s_h = int(h / factor)
    if logo_w > s_w or logo_h > s_h:
        logo_w = s_w
        logo_h = s_h

    logo = logo.resize((logo_w, logo_h), Image.ANTIALIAS)
    l_w = int((w - logo_w) / 2)
    l_h = int((h - logo_h) / 2)
    logo = logo.convert("RGBA")
    img.paste(logo, (l_w, l_h), logo)
    img.show()
    img.save(os.getcwd() + '/' + qrcodename + '.png', quality=100)


url = "http://2018.igem.org/Team:NEFU_China"
qrcodename = "igem_nefu2"
create_qrcode(url, qrcodename)