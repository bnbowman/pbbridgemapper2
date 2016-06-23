default:
	echo "nothing is specified, exiting"

pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'pbbridgemapper=='>/dev/null \
      && pip uninstall -y pbbridgemapper \
      || true
	@pip install --no-index \
          --install-option="--install-scripts=$(PREFIX)/bin" ./

