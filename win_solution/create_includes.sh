pushd ../src/;
ls ./*.cpp | awk '{print "    <ClCompile Include=\"..\\src\\" $1"\"\n      <Filter>Source Files</Filter>\n    </ClCompile>"}';
ls ./*.h | awk '{print "    <ClCompile Include=\"..\\src\\" $1"\"\n      <Filter>Header Files</Filter>\n    </ClCompile>"}';
ls ./*.cpp | awk '{print "    <ClCompile Include=\"..\\src\\" $1"\" />"}';
ls ./*.h | awk '{print "    <ClCompile Include=\"..\\src\\" $1"\" />"}';
