'''
===================================================================================================
Copyright (C) 2023 Steven Linfield

This file is part of the electrochemistry-simulations package. This package is free software: you 
can redistribute it and/or modify it under the terms of the GNU General Public License as published 
by the Free Software Foundation, either version 3 of the License, or (at your option) any later 
version. This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details. You should have received a copy of the GNU General 
Public License along with electrochemistry-simulations. If not, see https://www.gnu.org/licenses/
===================================================================================================

Package title:      electrochemistry-simulations
Repository:         https://github.com/MoonshinetheP/electrochemistry-simulations
Date of creation:   09/03/2023
Main author:        Steven Linfield (MoonshinetheP)
Collaborators:      None
Acknowledgements:   Oliver Rodriguez (oliverrdz)
    
Filename:           prediction.py

===================================================================================================
How to use this file:
    

===================================================================================================
'''

import tensorflow as tf
from keras.layers import Input, Dense
from keras.models import Model

num_input_features = 10

input_tensor = Input(shape=(2,))
hidden_layer = Dense(64, activation='relu')(input_tensor)
output_tensor = Dense(2)(hidden_layer)

model = Model(input_tensor, output_tensor)

model.compile(optimizer='adam', loss='mse')

num_epochs = 100
batch_size = 32
model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=num_epochs, batch_size=batch_size)
