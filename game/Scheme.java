package game;

import java.io.*;

// Loads a Bomberman scheme file.
public class Scheme
{
    // Special chars in file
    private static final char COMMENT_CHAR       = ';';
    private static final char DATA_SEPERATOR     = ',';

    private static final char CONTROL_CHAR       = '-';
    private static final char VERSION_CHAR       = 'V';
    private static final char NAME_CHAR          = 'N';
    private static final char BRICK_DENSITY_CHAR = 'B';
    private static final char DATA_ARRAY_CHAR    = 'R';
    private static final char PLAYER_START_CHAR  = 'S';
    private static final char POWER_UP_CHAR      = 'P';

    // Important scheme constants
    public static final char SOLID_CHAR = '#';
    public static final char BRICK_CHAR = ':';
    public static final char BLANK_CHAR = '.';

    // To read from file
    private BufferedReader in;

    // Read from file
    private int scheme_version;
    private String scheme_name;
    private int brick_density;
    // Standard sceme width = 15 and height = 11
    private char[][] scheme = new char[15][11];
    // Standard number of players = 10
    // players[player_number][0] = x
    // players[player_number][1] = y
    // players[player_number][2] = team
    private int[][] players = new int[10][3];
    // Standard number of powerups = 13
    // powerups[powerup_number][0] = born with (boolean)
    // powerups[powerup_number][1] = has override (boolean)
    // powerups[powerup_number][2] = override number
    // powerups[powerup_number][3] = forbidden (boolean)
    // NOT USED YET, UNKNOWN FUNCTION
    private int[][] powerups = new int[13][4];

    // Logging facility
    private PrintStream log;

    public Scheme(String file_path, String file_name, PrintStream log) throws Exception
    {
        this.log = log;
        this.in  = new BufferedReader(new FileReader(file_path + file_name));
        parseFile();
    }

    // Parse schemefile
    private void parseFile() throws Exception
    {
        String input;

        while((input = nextLine()) != null)
        {
            // Every line's first char is -, so ignore it
            if(input.charAt(0) != CONTROL_CHAR)
            {
                throw new Exception("Wrong control character in file (No -).");
            }

            char ctrl_char = input.charAt(1);
            input = input.substring(input.indexOf(DATA_SEPERATOR) + 1).trim();

            switch(ctrl_char)
            {
                case VERSION_CHAR:
                    String version;
                    writeLog(version = input);
                    this.scheme_version = Integer.valueOf(version);
                    break;
                case NAME_CHAR:
                    writeLog(this.scheme_name = input);
                    break;
                case BRICK_DENSITY_CHAR:
                    // Unknown what is does
                    String density;
                    writeLog(density = input);
                    this.brick_density = Integer.valueOf(density);
                    break;
                case DATA_ARRAY_CHAR:
                    // Get array index from data
                    String c_index = input.substring(0, input.indexOf(DATA_SEPERATOR)).trim();
                    // Remove index from input
                    input = input.substring(input.indexOf(DATA_SEPERATOR) + 1).trim();
                    int char_index = Integer.valueOf(c_index);

                    // Load data into array
                    for(int i = 0; i < input.length(); ++i)
                    {
                        this.scheme[i][char_index] = input.charAt(i);
                    }

                    break;
                case PLAYER_START_CHAR:
                    // Get player index from data
                    String p_index = input.substring(0, input.indexOf(DATA_SEPERATOR)).trim();
                    int player_index = Integer.valueOf(p_index);
                    // Remove index from input
                    input = input.substring(input.indexOf(DATA_SEPERATOR) + 1).trim();

                    // Get player x-pos from data
                    String x_pos = input.substring(0, input.indexOf(DATA_SEPERATOR)).trim();
                    int player_x_pos = Integer.valueOf(x_pos);
                    // Remove x-pos from input
                    input = input.substring(input.indexOf(DATA_SEPERATOR) + 1).trim();

                    // Get player y-pos from data
                    String y_pos = input.substring(0, input.indexOf(DATA_SEPERATOR)).trim();
                    int player_y_pos = Integer.valueOf(y_pos);
                    // Remove y-pos from input
                    input = input.substring(input.indexOf(DATA_SEPERATOR) + 1).trim();

                    // Get player team from data
                    int player_team = Integer.valueOf(input);

                    // Load data into this.players
                    this.players[player_index][0] = player_x_pos;
                    this.players[player_index][1] = player_y_pos;
                    this.players[player_index][2] = player_team;
                    break;
                case POWER_UP_CHAR:
                    // Get powerup index from data
                    String u_index = input.substring(0, input.indexOf(DATA_SEPERATOR)).trim();
                    int pu_index = Integer.valueOf(u_index);
                    // Remove index from input
                    input = input.substring(input.indexOf(DATA_SEPERATOR) + 1).trim();

                    // Get born with from data
                    String born_with = input.substring(0, input.indexOf(DATA_SEPERATOR)).trim();
                    int pu_born_with = Integer.valueOf(born_with);
                    // Remove born with from input
                    input = input.substring(input.indexOf(DATA_SEPERATOR) + 1).trim();

                    // Get has override from data
                    String has_override = input.substring(0, input.indexOf(DATA_SEPERATOR)).trim();
                    int pu_has_override = Integer.valueOf(has_override);
                    // Remove has override from input
                    input = input.substring(input.indexOf(DATA_SEPERATOR) + 1).trim();

                    // Get override number from data
                    String override_number = input.substring(0, input.indexOf(DATA_SEPERATOR)).trim();
                    int pu_override_number = Integer.valueOf(override_number);
                    // Remove override number from input
                    input = input.substring(input.indexOf(DATA_SEPERATOR) + 1).trim();

                    // Get forbidden from data
                    String is_forbidden;
                    int seperator_index = input.indexOf(DATA_SEPERATOR);

                    if (seperator_index < 0)
                    {
                        // No comment added
                        is_forbidden = input.trim();
                    }
                    else
                    {
                        // Comment added
                        is_forbidden = input.substring(0, input.indexOf(DATA_SEPERATOR)).trim();
                    }

                    int pu_is_forbidden = Integer.valueOf(is_forbidden);

                    // Load data into this.players
                    this.powerups[pu_index][0] = pu_born_with;
                    this.powerups[pu_index][1] = pu_has_override;
                    this.powerups[pu_index][2] = pu_override_number;
                    this.powerups[pu_index][3] = pu_is_forbidden;
                    break;
                default:
                    throw new Exception("Wrong control character in file (-" + ctrl_char + ").");
            }
        }
    }

    // Get trimmed next line from in
    private String nextLine() throws Exception
    {
        String input;

        // Ignore whitelines and comments
        while (((input = in.readLine()) != null) &&
                (input.trim().equals("") ||
                (input.trim().charAt(0) == COMMENT_CHAR)));

        return (input != null)?input.trim():null;
    }

    public int getSchemeVersion()
    {
        return this.scheme_version;
    }

    public int getBrickDensity()
    {
        return this.brick_density;
    }

    public String getSchemeName()
    {
        return this.scheme_name;
    }

    public char[][] getScheme()
    {
        return this.scheme;
    }

    public int[][] getPlayers()
    {
        return this.players;
    }

    public int[][] getPowerups()
    {
        return this.powerups;
    }

    // Logging facility
    private void writeLog(String s)
    {
        if(log != null)
        {
            log.println(s);
        }
    }

    // Testing facility
    public static void main(String[] args)
    {
        try
        {
            new Scheme("DATA/", "test.sch", System.out);
        }
        catch(Exception e)
        {
            e.printStackTrace();
            System.exit(0);
        }
    }
}