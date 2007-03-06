/**
 * A <code>InitCommand</code> contains data that is send by a server to a
 * newly connected client. The command contains the scheme name that is used
 * and the starting position of the client.
 */
public class InitCommand extends Command
{
    // Something
    public static final long serialVersionUID = 1;

    // The unique client id given by the server
    private int new_client_id;

    // The scheme that is played
    private String scheme_name;

    // The horizontal starting position of the client.
    private int x_pos;

    // The vertical starting position of the client.
    private int y_pos;

    /**
     * Creates a new <code>InitCommand</code> whith the given parameters.
     *
     * @param sender_id     The id of the sender of this command
     * @param new_client_id The unique id of the client to
     *                      which this command is send.
     * @param scheme_name   The scheme name that is played.
     * @param x_pos         The horizontal starting position of the client.
     * @param y_pos         The vertical starting position of the client.
     */
    public InitCommand(int sender_id, int new_client_id,
                       String scheme_name, int x_pos, int y_pos)
    {
        super(Config.CMD_INIT, sender_id);
        this.new_client_id = new_client_id;
        this.scheme_name = scheme_name;
        this.x_pos = x_pos;
        this.y_pos = y_pos;
    }

    /**
     * Returns the unique client id of given by the server.
     *
     * @return The unique client id given by the server.
     */
    public int getClientID()
    {
        return new_client_id;
    }

    /**
     * Returns the scheme name of this game.
     *
     * @return The scheme name.
     */
    public String getSchemeName()
    {
        return scheme_name;
    }

    /**
     * Returns the horizontal starting position of the client.
     *
     * @return The horizontal starting position of the client.
     */
    public int getXPos()
    {
        return this.x_pos;
    }

    /**
     * Returns the vertical starting position of the client.
     *
     * @return The vertical starting position of the client.
     */
    public int getYPos()
    {
        return this.y_pos;
    }
}