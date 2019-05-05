package ai.dragonfly.img.ops

import ai.dragonfly.color._
import ai.dragonfly.color.Color._
import ai.dragonfly.img.Img
import ai.dragonfly.math.stats.StreamingVectorStats
import ai.dragonfly.math.stats.kernel._
import ai.dragonfly.math.vector.{Vector2, Vector3, VectorN}

import scala.collection.mutable

package object ops {

  def randomizeRGB(img: Img): Img = {
    img pixels ( (x: Int, y: Int) => img.setARGB( x, y, Color.random().argb ) )
  }

  def randomizeLab(img: Img): Img = {
    img pixels ( (x: Int, y: Int) => img.setARGB( x, y, Color.randomFromLabSpace() ) )
  }

  def flipHorizontal(img: Img): Img = {
    val flipped: Img = new Img(img.width, img.height)
    val end:Int = img.width - 1
    flipped pixels ((x: Int, y: Int) => {
      flipped.setARGB( x, y, img.getARGB(end - x, y) )
    })
  }

  def flipVertical(img: Img): Img = {
    val flipped: Img = new Img(img.width, img.height)
    val end:Int = img.height - 1
    flipped pixels ((x: Int, y: Int) => {
      flipped.setARGB( x, y, img.getARGB(x, end - y) )
    })
  }

  def rotate90Degrees(img: Img, counterClockwise: Boolean = false): Img = {
    val rotated: Img = new Img(img.height, img.width)
    val rotationFunction = if (counterClockwise) {
      val endY:Int = img.width - 1; (x: Int, y: Int) => rotated.setARGB(x, y, img.getARGB(endY - y, x))
    } else {
      val endX = img.height - 1; (x: Int, y: Int) => rotated.setARGB(x, y, img.getARGB(y, endX - x))
    }
    rotated pixels rotationFunction
  }

  def rotate180Degrees (img: Img): Img = {
    val rotated: Img = new Img(img.width, img.height)
    val endX: Int = img.width - 1
    val endY: Int = img.height - 1
    rotated pixels ((x: Int, y: Int) => rotated.setARGB(x, y, img.getARGB(endX - x, endY - y)))
  }

  def overlay(
     bgImg: Img,
     fgImg: Img,
     bgX: Int, bgY: Int, fgX: Int, fgY: Int,
     width: Int, height: Int
   ): Img = {
    for (y <- 0 until height) {
      for (x <- 0 until width) {
        val fgc: RGBA = fgImg.getARGB(fgX + x, fgY + y)
        val bgc: RGBA = bgImg.getARGB(bgX + x, bgY + y)
        bgImg.setARGB(bgX + x, bgY + y, Color.alphaBlend(fgc, bgc).argb)
      }
    }
    bgImg
  }


  def epanechnikovBlurRGB(toBlur: Img, radius: Int): Img = kernelBlurRGB(toBlur, EpanechnikovKernel(radius))

  def uniformBlurRGB(toBlur: Img, radius: Int): Img = kernelBlurRGB(toBlur, UniformKernel(radius))

  def gaussianBlurRGB(toBlur: Img, radius: Int): Img = kernelBlurRGB(toBlur, GaussianKernel(radius))


  def kernelBlurRGB(toBlur: Img, kernel: Kernel): Img = {
    val dk = kernel.discretize

    val width:Int = toBlur.width
    val height:Int = toBlur.height

    val temp = new Img(width, height)
    val r: Int = Math.ceil(dk.radius).toInt

    val vectorStats = new StreamingVectorStats(3)

    // First Pass
    toBlur pixels ((x: Int, y: Int) => {
      for ( xi: Int <- Math.max(0, x - r) until Math.min(width, x + r + 1) ) {
        val dx = xi - x
        val c: RGBA = toBlur.getARGB(xi, y)
        vectorStats(new Vector3(c.red, c.green, c.blue), dk.weight(dx*dx))
      }
      val avg: Vector3 = vectorStats.average().asInstanceOf[Vector3]
      temp.setARGB(x, y, RGBA(avg.x.toInt, avg.y.toInt, avg.z.toInt).argb)
      vectorStats.reset()
    })

    toBlur pixels ((x: Int, y: Int) => {
      for ( yi <- Math.max(0, y - r) until Math.min(height, y + r + 1) ) {
        val dy = yi - y
        val c: RGBA = temp.getARGB(x, yi)
        vectorStats(new Vector3(c.red, c.green, c.blue), dk.weight(dy*dy))
      }
      val avg: Vector3 = vectorStats.average().asInstanceOf[Vector3]
      toBlur.setARGB(x, y, RGBA(avg.x.toInt, avg.y.toInt, avg.z.toInt).argb)
      vectorStats.reset()
    })
  }

  private def rgbaClamp(i: Double): Int = Math.max(Math.min(255, i), 0).toInt

  def unsharpenMaskRGB(img: Img, radius: Int, amount: Double, threshold: Int = 0): Img = {
    //sharpened = original + (original − blurred) × amount
    val t2 = threshold * threshold
    val blurred = gaussianBlurRGB(img.copy(), radius)

    img pixels ((x: Int, y: Int) => {
      val oc = img.getARGB(x, y)
      val bc = blurred.getARGB(x, y)

      var r = oc.red
      val rb = bc.red
      val rDif = r - rb

      if (rDif * rDif >= t2) r = rgbaClamp(r + (rDif * amount))

      var g = oc.green
      val gb = bc.green
      val gDif = g - gb
      if (gDif * gDif >= t2) g = rgbaClamp(g + (gDif * amount))

      var b = oc.blue
      val bb = bc.blue
      val bDif = b - bb
      if (bDif * bDif >= t2) b = rgbaClamp(b + (bDif * amount))

      img.setARGB(x, y, RGBA(r, g, b, oc.alpha))
    })
  }

  def unsharpenMaskLAB(img: Img, radius: Int, amount: Double, threshold: Int = 0): Img = {
    val t2 = threshold * threshold
    val blurred = gaussianBlurRGB(img.copy(), radius)

    img pixels ((x: Int, y: Int) => {
      val oc:LAB = RGBA(img.getARGB(x, y))
      val bc:LAB = RGBA(blurred.getARGB(x, y))

      var l1:Double = oc.L
      val l2:Double = bc.L
      val diff = l1 - l2

      if (diff * diff >= t2) l1 = l1 + (diff * amount)

      val c = SlowSlimLab(l1.toFloat, oc.a, oc.b).argb
      img.setARGB(x, y, RGBA(c.red, c.green, c.blue, oc.alpha))
    })
  }

  /*
   * Creates a difference matte from the two input images.
   *
   * @param img1: Img  The dimensions of this image determines the dimensions of the output image.
   * @param img2: Img  If this image does not have the same dimensions as img1, it is scaled to fit.
   * @return ImageBasics an image representing the difference matte.

   */

  def differenceMatte(img1: Img, img2: Img): Img = {

    val fitted = if ( img1.width != img2.width || img1.height != img2.height ) scale(img2, img1.width, img2.height)
    else img2

    val comparison: Img = new Img( fitted.width,  fitted.height )

    comparison pixels ((x: Int, y: Int) => {
      val c1: RGBA = img1.getARGB( x, y )
      val c2: RGBA = fitted.getARGB( x, y )

      // Compute the difference
      val dif = RGBA( Math.abs ( c1.red - c2.red ), Math.abs( c1.green - c2.green ), Math.abs ( c1.blue - c2.blue ) )
      comparison.setARGB(x, y, dif)
    })
  }

  // scale images.  Bilinear interpolation

  def scale(img: Img, newWidth: Int, newHeight: Int ): Img = {

    if (newWidth >= img.width && newHeight >= img.height) { // Bilinear interpolation to scale image up.
      val scaleX: Double = img.width / newWidth.toDouble
      val scaleY: Double = img.height / newHeight.toDouble

      val scaled: Img = new Img(newWidth, newHeight)

      scaled pixels ((u: Int, v: Int) => {
        val u1 = scaleX * u
        val v1 = scaleY * v

        val x1 = u1 - Math.floor(u * scaleX)
        val y1 = v1 - Math.floor(v * scaleY)

        val sU = u1.toInt
        val eU = Math.min(img.width - 1, u1 + 1).toInt
        val sV = v1.toInt
        val eV = Math.min(img.height - 1, v1 + 1).toInt

        val c00: RGBA = img.getARGB(sU, sV)
        val c01: RGBA = img.getARGB(sU, eV)
        val c10: RGBA = img.getARGB(eU, sV)
        val c11: RGBA = img.getARGB(eU, eV)

        val w1 = (1 - x1) * (1 - y1)
        val w2 = x1 * (1 - y1)
        val w3 = (1 - x1) * y1
        val w4 = x1 * y1

        val red: Int = (c00.red * w1 + c10.red * w2 + c01.red * w3 + c11.red * w4).toInt
        val green: Int = (c00.green * w1 + c10.green * w2 + c01.green * w3 + c11.green * w4).toInt
        val blue: Int = (c00.blue * w1 + c10.blue * w2 + c01.blue * w3 + c11.blue * w4).toInt

        scaled.setARGB(u, v, RGBA(red, green, blue))
      })
    } else if (newWidth <= img.width && newHeight <= img.height) {  // sampling to shrink image
      val scaleX = newWidth.toDouble / img.width
      val scaleY = newHeight.toDouble / img.height

      val statsImg: Array[Array[StreamingVectorStats]] = Array.fill(newWidth, newHeight){ new StreamingVectorStats(4)  }
      img pixels ((x: Int, y: Int) => {
        val c: RGBA = img.getARGB(x, y)
        statsImg((x * scaleX).toInt)((y * scaleY).toInt)(new VectorN(c.alpha, c.red, c.green, c.blue))
      })

      val scaled: Img = new Img(newWidth, newHeight)
      scaled pixels ((x: Int, y: Int) => {
        scaled.setARGB(x, y, {
          val v: Array[Double] = statsImg(x)(y).average().values
          RGBA(v(1).toInt, v(2).toInt, v(3).toInt, v(0).toInt).argb
        })
      })
    } else { // Shrink one dimension and grow the other

      val shrinkFirst = if (newWidth < img.width) scale(img, newWidth, img.height)
      else scale(img, img.width, newHeight)

      scale(shrinkFirst, newWidth, newHeight)
    }
  }

  def rotateDegrees(img: Img, angleDegrees: Double): Img = rotateRadians(img, angleDegrees * 0.01745329252)

  def rotateRadians(img: Img, angleRadians: Double): Img = {
    // Step 1, assess canvas size for resulting image:
    val midpoint1 = new Vector2(img.width / 2.0, img.height / 2.0)
    val corners: Array[Vector2] = Array(
      Vector2(0, 0).subtract(midpoint1),
      Vector2(0, img.height).subtract(midpoint1),
      Vector2(img.width, 0).subtract(midpoint1),
      Vector2(img.width, img.height).subtract(midpoint1)
    )

    var minX = Double.MaxValue
    var minY = Double.MaxValue
    var maxX = Double.MinValue
    var maxY = Double.MinValue

    for ( v <- corners ) {
      val rotated = v.rotate(angleRadians)
      minX = Math.min(rotated.x, minX)
      minY = Math.min(rotated.y, minY)
      maxX = Math.max(rotated.x, maxX)
      maxY = Math.max(rotated.y, maxY)
    }

    val rotated: Img = new Img(Math.sqrt(Math.pow(maxX - minX, 2)).toInt + 2, Math.sqrt(Math.pow(maxY - minY, 2)).toInt + 2)

    val midpoint2: Vector2 = new Vector2(rotated.width / 2.0, rotated.height / 2.0)

    // Sample pixels from img
    rotated pixels ((x: Int, y: Int) => {

      val x0 = x - midpoint2.x
      val y0 = y - midpoint2.y

      val cos = Math.cos( -angleRadians )
      val sin = Math.sin( -angleRadians )

      val u1 = (x0 * cos) + (y0 * -sin) + midpoint1.x
      val v1 = (x0 * sin) + (y0 * cos) + midpoint1.y

      if (u1 >= 0 && u1 < img.width && v1 >= 0 && v1 < img.height) {
        val x1 = u1 - Math.floor(u1)
        val y1 = v1 - Math.floor(v1)

        val sU = u1.toInt
        val eU = Math.min(img.width - 1, u1 + 1).toInt
        val sV = v1.toInt
        val eV = Math.min(img.height - 1, v1 + 1).toInt

        val c00: RGBA = img.getARGB(sU, sV)
        val c01: RGBA = img.getARGB(sU, eV)
        val c10: RGBA = img.getARGB(eU, sV)
        val c11: RGBA = img.getARGB(eU, eV)

        val w1 = (1 - x1) * (1 - y1)
        val w2 = x1 * (1 - y1)
        val w3 = (1 - x1) * y1
        val w4 = x1 * y1

        val red: Int = (c00.red * w1 + c10.red * w2 + c01.red * w3 + c11.red * w4).toInt
        val green: Int = (c00.green * w1 + c10.green * w2 + c01.green * w3 + c11.green * w4).toInt
        val blue: Int = (c00.blue * w1 + c10.blue * w2 + c01.blue * w3 + c11.blue * w4).toInt

        rotated.setARGB(x, y, RGBA(red, green, blue))
      }

    })
  }

  def grayscaleAverageRGB(img: Img): Img = {
    img.pixels((x: Int, y: Int) => {
      val c = img.getARGB(x, y)
      val avgIntensity = (c.red + c.green + c.blue) / 3
      img.setARGB(x, y, RGBA(avgIntensity, avgIntensity, avgIntensity, c.alpha))
    })
  }

  def grayscaleLABIntensity(img: Img): Img = {
    img.pixels((x: Int, y: Int) => {
      val c: RGBA = RGBA(img.getARGB(x, y))
      val lab: LAB = c
      val intensity: RGBA = SlowSlimLab(lab.L, 0f, 0f)
      img.setARGB(x, y, RGBA(intensity.red, intensity.green, intensity.blue, c.alpha))
    })
  }

  def negative(img: Img): Img = {
    img pixels ((x: Int, y: Int) => {
      val c = img.getARGB(x, y)
      img.setARGB(x, y, RGBA(255 - c.red, 255 - c.green, 255 - c.blue, c.alpha))
    })
  }

  def thresholdLab(img: Img, intensityRGB: Int): Img = {
    val labIntensity: LAB = RGBA(intensityRGB, intensityRGB, intensityRGB)
    val L: Double = labIntensity.L
    img pixels ((x: Int, y: Int) => {
      val c: RGBA = img.getARGB(x, y)
      val lab: LAB = c
      if (lab.L > L) img.setARGB(x, y, RGBA(255, 255, 255, c.alpha))
      else img.setARGB(x, y, RGBA(0, 0, 0, c.alpha))
    })
  }

  def thresholdRGB(img: Img, intensity: Int): Img = {
    img pixels ((x: Int, y: Int) => {
      val c: RGBA = img.getARGB(x, y)
      img.setARGB(x, y, RGBA(
        if (c.red > intensity) 255 else 0,
        if (c.green > intensity) 255 else 0,
        if (c.blue > intensity) 255 else 0,
        c.alpha
      ))
    })
  }

  def brightness(img: Img, b: Int): Img = {
    img pixels ((x: Int, y: Int) => {
      val c: RGBA = img.getARGB(x, y)
      img.setARGB(x, y, RGBA(
        Math.max(Math.min(255, c.red + b), 0),
        Math.max(Math.min(255, c.green + b), 0),
        Math.max(Math.min(255, c.blue + b), 0),
        c.alpha
      ))
    })
  }

  def contrast(img: Img, c: Double): Img = {
    val f = (259 * (c + 255)) / (255 * (259 - c))
    img pixels ((x: Int, y: Int) => {
      val c: RGBA = img.getARGB(x, y)
      img.setARGB(x, y, RGBA(
        Math.max(Math.min(255, (f * (c.red - 128) + 128).toInt), 0),
        Math.max(Math.min(255, (f * (c.green - 128) + 128).toInt), 0),
        Math.max(Math.min(255, (f * (c.blue - 128) + 128).toInt), 0),
        c.alpha
      ))
    })
  }

  def equalizeRGB(img: Img): Img = {
    val redHist = Array.fill[Double](256)(0.0)
    val greenHist = Array.fill[Double](256)(0.0)
    val blueHist = Array.fill[Double](256)(0.0)

    // compute histograms:
    img pixels ((x: Int, y: Int) => {
      val c = img.getARGB(x, y)
      redHist(c.red) = redHist(c.red) + 1.0
      greenHist(c.green) = greenHist(c.green) + 1.0
      blueHist(c.blue) = blueHist(c.blue) + 1.0
    })

    val pixelCount: Double = img.width * img.height

    var redCPD:Double = 0.0
    var greenCPD:Double = 0.0
    var blueCPD:Double = 0.0

    for (i <- 0 until 256) {
      redHist(i) = redCPD + (redHist(i) / pixelCount)
      redCPD = redHist(i)
      greenHist(i) = greenCPD + (greenHist(i) / pixelCount)
      greenCPD = greenHist(i)
      blueHist(i) = blueCPD + (blueHist(i) / pixelCount)
      blueCPD = blueHist(i)
    }

    img pixels ((x: Int, y: Int) => {
      val c = img.getARGB(x, y)
      img.setARGB(x, y, RGBA(
        Math.floor(redHist(c.red) * 255).toInt,
        Math.floor(greenHist(c.green) * 255).toInt,
        Math.floor(blueHist(c.blue) * 255).toInt,
        c.alpha
      ))
    })
  }

  def median(img: Img, radius: Int): Img = {
    val medianCut: Img = new Img(img.width, img.height)

    for (y <- 0 until img.height) {

      val r = new Histogram
      val g = new Histogram
      val b = new Histogram

      val startY = Math.max(0, y - radius)
      val endY = Math.min(medianCut.height - 1, y + radius)

      // initialize histograms with half window:
      for (y0 <- startY to endY) {
        for (x0 <- 0 to Math.min(medianCut.width - 1, radius)) {
          val c1 = img.getARGB(x0, y0)
          r(c1.red)
          g(c1.green)
          b(c1.blue)
        }
      }

      medianCut.setARGB(0, y, RGBA(r.median, g.median, b.median))

      for (x <- 1 until img.width) {

        val startX = Math.max(0, x - radius)
        val endX = Math.min(medianCut.width - 1, x + radius)

        // val count = (endY - startY) * (endX - startX)

        val oldColumn = startX - 1
        if (oldColumn >= 0) { // remove leftmost column
          for (y0 <- startY to endY) {
            val c = img.getARGB(oldColumn, y0)
            r.remove(c.red)
            g.remove(c.green)
            b.remove(c.blue)
          }
        }
        val newColumn = endX + 1
        if (newColumn < medianCut.width) {
          for (y0 <- startY to endY) {
            val c = img.getARGB(newColumn, y0)
            r(c.red)
            g(c.green)
            b(c.blue)
          }
        }

        medianCut.setARGB(x, y, RGBA(r.median, g.median, b.median))
      }
    }
    medianCut
  }

}

class Histogram {

  private val h = mutable.TreeMap[Int, Int]()
  private var total: Int = 0

  def apply(k: Int): Histogram = {
    h.get(k) match {
      case Some(count: Int) =>
        h.put(k, count + 1)
        total = total + 1
      case None => h.put(k, 1)
    }
    this
  }

  def remove(k: Int): Histogram = {
    h.get(k) match {
      case Some(count: Int) =>
        if (count <= 1) h.remove(k)
        else {
          h.put(k, count - 1)
          total = total - 1
        }
      case None =>
    }
    this
  }

  def median: Int = {
    var sum: Int = 0
    var candidate = -1
    val end = total / 2

    for ((k, count) <- h) {
      if (sum > end) return candidate
      else {
        candidate = k
        sum = sum + count
      }
    }
    candidate
  }

}
