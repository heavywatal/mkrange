for x in **.{h,cpp}; do sed -i 's/\r/\n/g' $x; done
for x in **.{h,cpp}; do nkf --in-place $x; done
for x in **.{h,cpp}; do sed -i 's/[ \t]*$//' $x; done
for x in **.{h,cpp}; do sed -Ei 's/<(\w+stream)\.h>/<\1>/g' $x; done
for x in **.{h,cpp}; do sed -Ei 's/"(list).h"/<\1>/g' $x; done
for x in **.{h,cpp}; do sed -Ei 's/<(\w+)\.h>/<c\1>/g' $x; done
for x in **.{h,cpp}; do sed -Ei 's/<c(unistd)>/<\1.h>/g' $x; done
for x in **.{h,cpp}; do sed -i 's/list</std::list</g' $x; done
for x in **.{h,cpp}; do sed -i 's/using namespace std.*$//' $x; done
for x in **.{h,cpp}; do sed -Ei 's/(cout|endl|ifstream)/std::\1/g' $x; done
