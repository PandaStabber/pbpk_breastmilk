<!-- Terminal steps to use the model from scratch! -->
# Step 1: install xcode tools 
xcode-select --install

# Step 2: install homebrew package manager
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# Step 3: install pip package manager
sudo easy_install pip

# Step 4: install scipy
sudo easy_install scipy

# Step 5: install virtualenv
sudo pip install virtualenv

# Step 6: create virtual environment
cd /eth_github_job_env/ <!--/ move to the director - whatever it's called!/--> 
mkdir eth_github_job_env
virtualenv eth_github_job_env

# Step 7: activate virtual environment
cd eth_github_job_env
source bin/activate <!--/ you're now in the virtual environment!/-->




