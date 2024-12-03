# Main OS Lib Import
import os

# Main Class Start
class Clear_String:

    def __init__(self):

        Dir_Current = 'os.getcwd()'

        Want_To_Remove_String = 'NNNNNNNN'

        List_ParseForFile = self.Func__GetFileList(Dir_Current); del(Dir_Current);

        # Remove String...
        self.Func__Remove(List_ParseForFile, Want_To_Remove_String)

    def Func__GetFileList(self, Dir_Current):
        List_File = os.listdir(Dir_Current)
        List_Result = []

        for x in List_File:

            if 'fastp.dj.fastq' in x:

                List_Result.append(x);

            del(x)

        del(Dir_Current); del(List_File);

        return(List_Result)

    def Func__Remove(self, List_ParseForFile, Want_To_Remove_String):

        for FileName in List_ParseForFile:

            FR = open(FileName, 'r', encoding='utf-8'); Lines = FR.read(); FR.close(); del(FR);

            ModifiedLines = Lines.replace(Want_To_Remove_String, '');

            FW = open(FileName, 'w', encoding='utf-8'); FW.write(ModifiedLines); FW.close(); del(FW);

            print(f"{FileName} Complete!!")

            del(FileName);

Clear_String()
