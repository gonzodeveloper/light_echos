from keras.layers import *
from keras.models import *


def build_model(shape):
    model = Sequential()
    model.add(Conv2D(filters=64, kernel_size=3, activation='relu',
                     kernel_initializer='he_normal', padding='same'))