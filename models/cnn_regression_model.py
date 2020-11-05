from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D, MaxPooling2D
from tensorflow.keras.layers import Activation, Dropout, Flatten, Dense, BatchNormalization

input_shape = (40,40,1)

def create_model():
    model = Sequential()
    model.add(Conv2D(32, kernel_size=(3, 3),
					 activation='relu',
                     padding='same',
                     input_shape=input_shape))
    model.add(BatchNormalization())
    model.add(Conv2D(32, kernel_size=(3, 3),
                         activation='relu',
                         padding='same',
                         name='block2_conv2'))
    model.add(BatchNormalization())
    model.add(MaxPooling2D(pool_size=(2, 2),
                         strides=(2, 2),
                         name='block2_pool'))
    model.add(Dropout(0.25))

    model.add(Conv2D(64, kernel_size=(3, 3),
                         activation='relu',
                         padding='same',
                         name='block3_conv1'))
    model.add(BatchNormalization())
    model.add(Conv2D(64, kernel_size=(3, 3),
                         activation='relu',
                         padding='same',
                         name='block3_conv2'))
    model.add(BatchNormalization())
    model.add(Conv2D(64, kernel_size=(3, 3),
                         activation='relu',
                         padding='same',
                         name='block3_conv3'))
    model.add(BatchNormalization())
    model.add(MaxPooling2D(pool_size=(2, 2),
                         strides=(2, 2),
                         name='block3_pool'))
    model.add(Dropout(0.25))

    model.add(Flatten())
    model.add(Dense(128, activation='relu'))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(32, activation='relu'))
    model.add(Dense(1, activation='relu'))
    model.add(Activation('linear'))

    return model






