'''
Created on 14.02.2015

@author: David
'''
from OptimizerApp import OptimizerApp
from OptimizerMainWindow import OptimizerMainWindow
import SeqUtils
import os, sys
import PathUtils

if __name__ == '__main__':
    
#     print(SeqUtils.getRemainderSuffix("AAAACCA"))
    configFile = os.path.join(PathUtils.getCwd(), "config.ini")
    print(configFile)
    myOptimizer = OptimizerApp(configFile)
    myOptimizerGUI = OptimizerMainWindow(myOptimizer)
#     gui = Tk()
#     gui.mainloop()
#     myOptimizer.testPrint()