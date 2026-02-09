#This is NOT an executable script. Each block that starts with comments should be run independently to account for potential errors. 

#Mount volume –– this step shouldn't be run in bulk with the rest. 
lsblk #assume volume is /dev/vdc
mkdir -p vol && sudo mount /dev/vdc vol 
sudo chmod 777 vol

#Directory where executables / scripts will be placed
mkdir -p ~/bin

#After all installations are complete and scripts are added, 
#Run from your machine:
rsync -avzP binderScore ubuntu@185.216.21.176:~/bin/

#Run from host machine:
ln -sf ~/bin/binderScore/binderScore.sh ~/bin/binderScore.sh
chmod +x ~/bin/*
chmod +x ~/bin/*/*
chmod +x ~/bin/*/*/*

echo 'export PATH="~/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

#conda install
wget -q "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-$(uname)-$(uname -m).sh
$HOME/miniforge3/bin/conda init bash
source ~/.bashrc
conda config --set auto_activate false
source ~/.bashrc

#colabfold install
curl -fsSL https://pixi.sh/install.sh | sh
source ~/.bashrc
git clone https://github.com/yoshitakamo/localcolabfold.git
cd localcolabfold
pixi install && pixi run setup
echo 'export PATH="/home/ubuntu/localcolabfold/.pixi/envs/default/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
conda create -n colabfold -c nvidia cuda-nvcc=12.4 -y
conda deactivate

rm -f ~/*.sh

#vmtouch install
git clone https://github.com/hoytech/vmtouch.git
cd vmtouch && make
sudo make install
cd ~

#jq install
sudo apt-get update
sudo apt install -y jq

#chai-1 install
conda create -n chai python==3.10 -y
conda activate chai
pip install --upgrade pip
pip install chai_lab==0.6.1
conda deactivate

#Boltz-2 install
conda create -n boltz python=3.11 -y
conda activate boltz
pip install boltz[cuda] -U
conda deactivate

#PyRosetta, dssp, DAlphaBall install
conda create -n pyrosetta python=3.11 -y
conda activate pyrosetta
pip install pyrosetta-installer
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(distributed=True,serialization=True)'
pip install pyrosetta-distributed\

conda install dssp -y
pip install biopython

git clone https://github.com/outpace-bio/DAlphaBall.git
mkdir DAlphaBall_build
./DAlphaBall/dalphaball_docker_build.sh ~/DAlphaBall_build
conda install dalphaball -c ~/DAlphaBall_build -y

conda deactivate



