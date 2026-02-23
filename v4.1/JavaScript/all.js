function switchTheme() {
  myLogo = document.getElementById("cassiopeeLogo")
  myBtn = document.getElementById("fancyButton")

  if (myBtn.textContent == "Add some yellow"){
    myLogo.style.backgroundColor = '#F5F5A0'
    myBtn.textContent = "Switch back"
  } 
  else {
    myLogo.style.backgroundColor = '#f0fdff'
    myBtn.textContent = "Add some yellow"
  } 
}

//  IMG-BOX - credit to: krittanon-w (https://www.cssscript.com/fullscreen-image-viewer-lightbox)

var bg_color_img_box = 'rgba(0,0,0,0.9)'
var allow_hide_scroll_img_box = 'yes'
var use_fade_inout_img_box = 'yes'
var speed_img_box = 0.08
var z_index_dv_img_box = 999
var vopa_img_box, idpopup_img_box

window.onload = function() {
    var crtdv_img_box = document.createElement('div')
    crtdv_img_box.id = 'img_box'
    document.getElementsByTagName('body')[0].appendChild(crtdv_img_box)
    idpopup_img_box = document.getElementById("img_box")
    idpopup_img_box.style.top = 0
    idpopup_img_box.style.left = 0
    idpopup_img_box.style.opacity = 0
    idpopup_img_box.style.width = '100%'
    idpopup_img_box.style.height = '100%'
    idpopup_img_box.style.display = 'none'
    idpopup_img_box.style.position = 'fixed'
    idpopup_img_box.style.cursor = 'pointer'
    idpopup_img_box.style.textAlign = 'center'
    idpopup_img_box.style.zIndex = z_index_dv_img_box
    idpopup_img_box.style.backgroundColor = bg_color_img_box
}

function img_box(self) {
    var namepic_img_box = typeof self === 'string' ? self : self.src
    vopa_img_box = 0
    var hwin_img_box = window.innerHeight
    var wwin_img_box = window.innerWidth
    var himg_img_box, wimg_img_box, padtop_img_box, idfadein_img_box
    var img_img_box = new Image()
    img_img_box.src = namepic_img_box
    img_img_box.onload = function() {
        himg_img_box = img_img_box.height
        wimg_img_box = img_img_box.width
        idpopup_img_box.innerHTML = '<img src=' + namepic_img_box + '>'

        idpopup_img_box.getElementsByTagName('img')[0].style.backgroundColor = '#ffffff';

        if (wimg_img_box > wwin_img_box) {
            idpopup_img_box.getElementsByTagName('img')[0].style.width = '80%'
            idpopup_img_box.getElementsByTagName('img')[0].style.height = 'auto'
            himg_img_box = (himg_img_box/wimg_img_box) * wwin_img_box * 80 / 100
        }
        else if (himg_img_box > hwin_img_box) {
            idpopup_img_box.getElementsByTagName('img')[0].style.height = '80%'
            idpopup_img_box.getElementsByTagName('img')[0].style.width  = 'auto'
            himg_img_box = hwin_img_box * 80 / 100
        }

        padtop_img_box = (hwin_img_box - himg_img_box) / 2
        idpopup_img_box.style.paddingTop = padtop_img_box + 'px'

        if (allow_hide_scroll_img_box == 'yes') {
            document.body.style.overflow = 'hidden'
        }
        idpopup_img_box.style.display = 'block'
    }

    if (use_fade_inout_img_box == 'yes') {
        idfadein_img_box = setInterval(function() {
            if (vopa_img_box <= 1.1) {
                idpopup_img_box.style.opacity = vopa_img_box
                vopa_img_box += speed_img_box
            }
            else {
                idpopup_img_box.style.opacity = 1
                clearInterval(idfadein_img_box)
            }
        }, 10)
    }
    else {
        idpopup_img_box.style.opacity = 1
    }

    idpopup_img_box.onclick = function() {
        if (use_fade_inout_img_box == 'yes') {
            var idfadeout_img_box = setInterval(function() {
                if (vopa_img_box >= 0) {
                    idpopup_img_box.style.opacity = vopa_img_box
                    vopa_img_box -= speed_img_box
                } else {
                    idpopup_img_box.style.opacity = 0
                    clearInterval(idfadeout_img_box)
                    idpopup_img_box.style.display = 'none'
                    idpopup_img_box.innerHTML = ''
                    document.body.style.overflow = 'visible'
                    vopa_img_box = 0
                }
            }, 10)
        }
        else {
            idpopup_img_box.style.opacity = 0
            idpopup_img_box.style.display = 'none'
            idpopup_img_box.innerHTML = ''
            document.body.style.overflow = 'visible'
        }
    }
}