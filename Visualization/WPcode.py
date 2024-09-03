# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 13:25:45 2024

@author: Lenovo
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 16:21:28 2024

@author: Lenovo
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 17:55:40 2024

@author: Lenovo
"""
import wordcloud
from wordcloud import WordCloud, ImageColorGenerator
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import imageio
import numpy as np

colormap = cm.get_cmap("Blues")
# 再把这个放到WordCloud()中的colormap参数值就行了

mk = imageio.imread("grey_gradiant2.png")

x, y = np.ogrid[:600, :600]
mask = (x - 300) ** 2 + (y - 300) ** 2 > 260 ** 2
mask = 255 * mask.astype(int)

from wordcloud import get_single_color_func


class SimpleGroupedColorFunc(object):
    def __init__(self, color_to_words, default_color):
        self.word_to_color = {word: color
                              for (color, words) in color_to_words.items()
                              for word in words}

        self.default_color = default_color

    def __call__(self, word, **kwargs):
        return self.word_to_color.get(word, self.default_color)

def create_word_cloud():
    frequencies = {}
    for line in open("wdata2.txt"):
        arr = line.split(",")
        frequencies[arr[0]] = float(arr[1])

    wc = WordCloud(
        prefer_horizontal = 1,
        mask = mask,
        font_path = "calibri.ttf",
        background_color = "white",
        contour_color = "yellow",
        colormap = colormap,
        max_words = 50,
        width = 600,
        height = 400)
    word_cloud = wc.generate_from_frequencies(frequencies)

    #color_to_words = {
     #   '#282a62': [word for word in frequencies.keys() if 'fear' in word or  'emot' in word
     #           or 'depress' in word or 'recall' in word or 'DMN' in word 
      #          or 'affect' in word or 'nega' in word or 'avoid' in word
     #           or 'memor' in word or 'detail' in word or 'hippo' in word
      #          or 'therapy' in word or 'BD' in word or 'MDD' in word
      #          or 'painful' in word or 'unpleas' in word
      #          or 'autobio' in word or 'bd' in word or 'bipolar' in word] 
  # }
#    color_to_words = {
#       '#c85d4d': [word for word in frequencies.keys() if 'risk' in word or  'negative' in word
 #              or 'drug' in word  or 'fear extinction' in word 
 #               or 'social interaction' in word or 'affective' in word or 'BD' in word
  #             or 'competition' in word or 'motor' in word or 'goal' in word
  #             or 'social interaction' in word or 'BD' in word
   #            or 'attention' in word or 'impul' in word or 'esteem' in word
   #            or 'PBD' in word or 'conflict' in word
   #            or 'decisiom' in word or 'interaction' in word
   #            or 'emotion' in word or 'decision' in word] 
   #}


#    default_color = 'Blues_r'
#    grouped_color_func = SimpleGroupedColorFunc(color_to_words, default_color)
    #grouped_color_func = SimpleGroupedColorFunc(default_color)

    #wc.recolor(color_func=grouped_color_func)

    word_cloud.to_file("conj.jpg")
    plt.imshow(word_cloud)
    plt.axis("off")
    plt.show()

p = create_word_cloud()