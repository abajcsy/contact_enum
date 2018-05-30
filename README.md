# Running on Docker

```bash
xhost +local:root; docker run -i --rm -e DISPLAY -e QT_X11_NO_MITSHM=1 -v /tmp/.X11-unix:/tmp/.X11-unix -v /path/to/contact_enum:/contact_enum --privileged -t eratner/drake_ellis; xhost -local:root
```
where /path/to/contact_enum should be replaced with the path to the contact_enum repository. 

```bash
cd /contact_enum
```

```bash
python src/run.py
```