#!/bin/bash

echo $(which conda)
conda init bash
conda activate python
echo $(which R)
conda deactivate
