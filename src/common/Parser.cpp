/*
 *================================================================================
 *
 *    Copyright (c) 2016 Vortex Co.,Ltd.
 *    Unpublished - All rights reserved
 *
 *================================================================================
 *    File description:
 *    Interface class for processor data
 *
 *================================================================================
 *    Date            Name                    Description of Change
 *    04-Oct-2016     Jiamin Xu               Creation
 *================================================================================
 */



using namespace std;


namespace ARIES
{

    Parser::Parser(string val_configFile):
            d_configFile(val_configFile)
    {
    }


    Parser::~Parser()
    {
    }
    
    void StringToUpperCase(string& str)
    {
        std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    }


    string StringToUpperCase(const string& str)
    {
        string upp_str(str);
        std::transform(upp_str.begin(), upp_str.end(), upp_str.begin(), ::toupper);
        return upp_str;
    }

    bool Parser::DoParse()
    {
        string textLine, optionName;
        ifstream caseFile;
        vector<string> optionValue;
        string errorString;
        int errCount = 0;      // How many errors have we found in the config file
        int maxErrCount = 30;  // Maximum number of errors to print before stopping

        int rank = MASTER_NODE;
  
#ifdef ARIES_HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        
        /*
         *  open the config file
         */
        caseFile.open(d_configFile, ios::in);

        if (caseFile.fail())
        {
            if (rank == MASTER_NODE)
                cout << endl << "The configuration file (.cfg) is missing!!" << endl << endl;
            exit(EXIT_FAILURE);
        }

        map<string, bool> included_options;

        /*--- Parse the configuration file and set the options ---*/
        while (getline (caseFile, textLine))
        {
            if (errCount >= maxErrCount)
            {
                errorString.append("too many errors. Stopping parse");
                cout << errorString << endl;
                throw(1);
            }
    
            if (TokenizeString(textLine, optionName, optionValue))
            {
                /*--- See if it's a python option ---*/
                if (option_map.find(optionName) == option_map.end())
                {
                    string newString;
                    newString.append(optionName);
                    newString.append(": invalid option name");
                    newString.append(". Check current SU2 options in config_template.cfg.");
                    newString.append("\n");
                    errorString.append(newString);
                    errCount++;
                    continue;
                }

                /*--- Option exists, check if the option has already been in the config file ---*/     
                if (included_options.find(optionName) != included_options.end())
                {
                    string newString;
                    newString.append(optionName);
                    newString.append(": option appears twice");
                    newString.append("\n");
                    errorString.append(newString);
                    errCount++;
                    continue;
                }

                /*--- New found option. Add it to the map, and delete from all options ---*/
                included_options.insert(pair<string, bool>(optionName, true));
                all_options.erase(optionName);

                /*--- Set the value and check error ---*/
                string out = option_map[optionName]->SetValue(optionValue);
                if (out.compare("") != 0)
                {
                    errorString.append(out);
                    errorString.append("\n");
                    errCount++;
                }
            }
        }

        /*--- See if there were any errors parsing the config file ---*/
        if (errorString.size() != 0)
        {
            if (rank == MASTER_NODE)
                cout << errorString << endl;
            exit(EXIT_FAILURE);
        }

        /*--- Set the default values for all of the options that weren't set ---*/
        for (map<string, bool>::iterator iter = all_options.begin(); iter != all_options.end(); ++iter)
        {
            option_map[iter->first]->SetDefault();
        }

        caseFile.close();
    }

    bool Parser::TokenizeString(string& str, string& optionName, vector<string>& optionValue)
    {
        const string delimiters(" ()[]{}:,\t\n\v\f\r");
        // check for comments or empty string
        string::size_type pos, lastPos;
        pos = str.find_first_of("%");
        if ( (str.length() == 0) || (pos == 0) )
        {
            // str is empty or a comment line, so no option here
            return false;
        }
        if (pos != string::npos)
        {
            // remove comment at end if necessary
            str.erase(pos);
        }

        // look for line composed on only delimiters (usually whitespace)
        pos = str.find_first_not_of(delimiters);
        if (pos == string::npos)
        {
            return false;
        }

        // find the equals sign and split string
        string name_part, value_part;
        pos = str.find("=");
        if (pos == string::npos)
        {
            cerr << "Error in TokenizeString(): "
                 << "line in the configuration file with no \"=\" sign."
                 << endl;
            cout << "Look for: " << str << endl;
            cout << "str.length() = " << str.length() << endl;
            throw(-1);
        }
        name_part = str.substr(0, pos);
        value_part = str.substr(pos+1, string::npos);
        //cout << "name_part  = |" << name_part  << "|" << endl;
        //cout << "value_part = |" << value_part << "|" << endl;

        // the first_part should consist of one string with no interior delimiters
        lastPos = name_part.find_first_not_of(delimiters, 0);
        pos = name_part.find_first_of(delimiters, lastPos);
        if ( (name_part.length() == 0) || (lastPos == string::npos) )
        {
            cerr << "Error in CConfig::TokenizeString(): "
                 << "line in the configuration file with no name before the \"=\" sign."
                 << endl;
            throw(-1);
        }
        if (pos == string::npos) pos = name_part.length();
        optionName = name_part.substr(lastPos, pos - lastPos);
        lastPos = name_part.find_first_not_of(delimiters, pos);
        if (lastPos != string::npos)
        {
            cerr << "Error in TokenizeString(): "
                 << "two or more options before an \"=\" sign in the configuration file."
                 << endl;
            throw(-1);
        }
        StringToUpperCase(optionName);

        //cout << "optionName = |" << optionName << "|" << endl;
        //cout << "pos = " << pos << ": lastPos = " << lastPos << endl;

        // now fill the option value vector
        optionValue.clear();
        lastPos = value_part.find_first_not_of(delimiters, 0);
        pos = value_part.find_first_of(delimiters, lastPos);
        while (string::npos != pos || string::npos != lastPos)
        {
            // add token to the vector<string>
            optionValue.push_back(value_part.substr(lastPos, pos - lastPos));
            // skip delimiters
            lastPos = value_part.find_first_not_of(delimiters, pos);
            // find next "non-delimiter"
            pos = value_part.find_first_of(delimiters, lastPos);
        }
        if (optionValue.size() == 0)
        {
            cerr << "Error in TokenizeString(): "
                 << "option " << optionName << " in configuration file with no value assigned."
                 << endl;
            throw(-1);
        }

#if 0
        cout << "option value(s) = ";
        for (unsigned int i = 0; i < optionValue.size(); i++)
            cout << optionValue[i] << " ";
        cout << endl;
#endif

        // look for ';' DV delimiters attached to values
        vector<string>::iterator it;
        it = optionValue.begin();
        while (it != optionValue.end())
        {
            if (it->compare(";") == 0)
            {
                it++;
                continue;
            }

            pos = it->find(';');
            if (pos != string::npos)
            {
                string before_semi = it->substr(0, pos);
                string after_semi= it->substr(pos+1, string::npos);
                if (before_semi.empty())
                {
                    *it = ";";
                    it++;
                    optionValue.insert(it, after_semi);
                }
                else
                {
                    *it = before_semi;
                    it++;
                    vector<string> to_insert;
                    to_insert.push_back(";");
                    if (!after_semi.empty())
                        to_insert.push_back(after_semi);
                    optionValue.insert(it, to_insert.begin(), to_insert.end());
                }
                it = optionValue.begin(); // go back to beginning; not efficient
                continue;
            }
            else
            {
                it++;
            }
        }
#if 0
        cout << "option value(s) = ";
        for (unsigned int i = 0; i < optionValue.size(); i++)
            cout << optionValue[i] << " ";
        cout << endl;
#endif
        // remove any consecutive ";"
        it = optionValue.begin();
        bool semi_at_prev = false;
        while (it != optionValue.end())
        {
            if (semi_at_prev)
            {
                if (it->compare(";") == 0)
                {
                    optionValue.erase(it);
                    it = optionValue.begin();
                    semi_at_prev = false;
                    continue;
                }
            }
            if (it->compare(";") == 0)
            {
                semi_at_prev = true;
            }
            else
            {
                semi_at_prev = false;
            }
            it++;
        }

#if 0
        cout << "option value(s) = ";
        for (unsigned int i = 0; i < optionValue.size(); i++)
            cout << optionValue[i] << " ";
        cout << endl;
#endif
        return true;
    }

}






























