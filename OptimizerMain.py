'''
Created on 14.02.2015

@author: David
'''
from OptimizerApp import OptimizerApp
from OptimizerMainWindow import OptimizerMainWindow
import SeqUtils


if __name__ == '__main__':
    
#     print(SeqUtils.getRemainderSuffix("AAAACCA"))
    
    myOptimizer = OptimizerApp("config.ini")
    myOptimizerGUI = OptimizerMainWindow(myOptimizer)
#     gui = Tk()
#     gui.mainloop()
#     myOptimizer.testPrint()